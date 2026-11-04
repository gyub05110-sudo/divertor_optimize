"""Check the saved display geometry, with no tracing or thermal evaluation."""
from export_data import HERE, read, sha, save
import numpy as np
from shapely.geometry import LineString, Point


def main():
    data = read(HERE/'dist/data.json')
    legs = read(HERE/'dist/legs.json')
    audit = read(HERE/'legs-audit.json')
    assert sha(HERE/'dist/data.json') == legs['data_sha256']
    assert sha(HERE/'dist/legs.json') == audit['display_leg_file_sha256']
    assert all(sha(path) == expected for path, expected in audit['source_sha256'].items())
    walls = {'baseline': data['baseline']['wall']}
    for series in data['series']:
        for row in series['stages'] + ([series['final']] if series['final'] else []):
            walls[row['id']] = row['wall']
    assert set(walls) == set(legs['cases'])
    count = 0
    for identity, pair in legs['cases'].items():
        wall = LineString(walls[identity])
        assert set(pair) == {'inner', 'outer'}
        for leg in pair.values():
            assert leg['resolved']
            points = np.asarray(leg['points'])
            assert points.ndim == 2 and points.shape[1] == 2 and len(points) >= 2
            assert np.all(np.isfinite(points))
            np.testing.assert_allclose(points[0], legs['xpoint_RZ_m'], atol=1e-12, rtol=0)
            np.testing.assert_allclose(points[-1], leg['strike_RZ_m'], atol=1e-12, rtol=0)
            assert wall.distance(Point(points[-1])) < 1e-9
            assert 0 < LineString(points).length <= leg['poloidal_length_m'] + 1e-12
            count += 1
    baseline = legs['cases']['baseline']['outer']
    omp50 = next(s['final'] for s in data['series'] if s['source']=='OMP' and s['rmin']==50 and s['run']=='tangent')
    selected = legs['cases'][omp50['id']]['outer']
    # The lowered selected plate must extend the branch beyond the original strike.
    assert selected['strike_RZ_m'][1] < baseline['strike_RZ_m'][1] - .02
    assert selected['poloidal_length_m'] > baseline['poloidal_length_m']
    report = {'checks':'passed', 'wall_cases':len(walls), 'leg_paths':count,
              'all_candidates_covered':True, 'all_endpoints_on_displayed_wall':True,
              'protected_input_hashes_unchanged':True, 'new_thermal_solves':0,
              'OMP50_outer_strike_baseline_m':baseline['strike_RZ_m'],
              'OMP50_outer_strike_selected_m':selected['strike_RZ_m'],
              'OMP50_outer_poloidal_length_baseline_m':baseline['poloidal_length_m'],
              'OMP50_outer_poloidal_length_selected_m':selected['poloidal_length_m']}
    save(HERE/'legs-verification.json', report)
    print(report)


if __name__ == '__main__':
    main()
