"""Opt-in endpoint broadening; the historical direct evaluator is untouched."""
from dataclasses import asdict
import numpy as np
from shapely.geometry import Point
from .broadening import BroadeningConfig, PathLengths, build_deposition
from .frozen_mapping import Hits

MODEL_VERSION = 'connection-deposition-surrogate.v1'


class BroadenedEvaluator:
    def __init__(self, direct, metadata, config: BroadeningConfig, lengths=None):
        self.direct, self.metadata, self.config = direct, dict(metadata), config
        self.path_lengths = (lengths if lengths is not None else PathLengths(direct.field, direct.paths)) if config.sigma_ref_m else None
        self.baseline = self.evaluate(np.zeros(direct.dimension))
        self.baseline_peaks = self.baseline.get('peak_q')
        self.baseline_admissible = bool(self.baseline['feasible'] and self.baseline_peaks is not None
                                        and np.all(np.isfinite(self.baseline_peaks))
                                        and np.all(self.baseline_peaks > 0))

    def __getattr__(self, name):
        return getattr(self.direct, name)

    def evaluate(self, coefficients):
        result = self.direct.evaluate(coefficients)
        if self.config.sigma_ref_m == 0:
            return result
        r = dict(result)
        r['model_version'] = MODEL_VERSION
        r['broadening_config'] = asdict(self.config)
        if not r['geometry_valid']:
            return r
        # Keep direct geometry/incidence evidence, but never mix direct peaks
        # or the zero-power separatrix heat witness into broadened heat peaks.
        r['direct_boundary_witness'] = r.pop('boundary_witness', None)
        r['direct_maps'] = r.pop('maps')
        r.pop('peak_q', None)
        hits = Hits(r['direct_maps'][0].segment, r['direct_maps'][0].points,
                    np.array([path.project(Point(point)) if seg >= 0 else np.inf
                              for path, point, seg in zip(self.direct.paths,
                                  r['direct_maps'][0].points, r['direct_maps'][0].segment)]))
        try:
            dep = build_deposition(
                self.direct.field, self.direct.paths, self.path_lengths, hits, r['wall'],
                np.column_stack([s.weights for s in self.direct.sources]),
                self.metadata['source_RZ_m'], -1 if self.metadata['source'] == 'IMP' else 1,
                self.direct.lcfs, self.direct.vessel_solids, self.config)
        except (ValueError, FloatingPointError) as error:
            r.update(feasible=False, penalty=1e6, model_errors=[str(error)])
            return r
        r['deposition'] = dep
        r['peak_q'] = np.max(dep['q'], axis=0)
        ledger = dep['ledger']
        r['missing'] = float(np.max(ledger['unresolved']))
        r['power_balance_error'] = float(np.max(ledger['balance_error']))
        finite = bool(np.all(np.isfinite(dep['q'])) and np.all(dep['q'] >= 0)
                      and np.all(np.isfinite(r['peak_q'])))
        receiving_angles = {}
        angle_ok = True
        angle_penalty = 0.
        incidence = self.direct.field.incidence(dep['points'], dep['normals'])
        for surface in self.direct.surfaces:
            name, floor = surface['name'], surface['angle_floor_deg']
            ids = r['surface_segments'][name]
            # A grazing endpoint has zero Jacobian/q even on a powered cell.
            # Use integrated receiving power, never q, for the angle policy.
            powered = np.any(dep['cell_power'][dep['sample_cell']] > 0, axis=1)
            mask = np.isin(dep['sample_parent'], ids) & powered
            minimum = float(np.min(incidence[mask])) if np.any(mask) else None
            passed = floor is None or minimum is None or (np.isfinite(minimum) and minimum >= floor)
            receiving_angles[name] = dict(minimum_angle_deg=minimum, floor_deg=floor,
                                         receiving_sample_count=int(mask.sum()), passed=bool(passed))
            if not passed:
                angle_ok = False
                angle_penalty += 100 * max(0., floor - minimum) if np.isfinite(minimum) else 1000.
        r['receiving_angles'] = receiving_angles
        r['surface_power'] = {
            name: np.sum(dep['cell_power'][np.isin(dep['cell_parent'], ids)], axis=0).tolist()
            for name, ids in r['surface_segments'].items()}
        illuminated = any(any(value > 0 for value in row) for row in r['surface_power'].values())
        r['feasible'] = bool(result['feasible'] and dep['valid'] and finite and angle_ok
                             and illuminated and r['missing'] <= self.config.missing_tolerance
                             and r['power_balance_error'] <= 1e-10)
        r['penalty'] = (result['penalty'] + angle_penalty
                        + 1000 * max(0., r['missing'] - self.config.missing_tolerance)
                        + (0. if dep['valid'] and finite else 1000.)
                        + (0. if illuminated else 1000.))
        r['model_errors'] = list(dep.get('reasons', []))
        return r
