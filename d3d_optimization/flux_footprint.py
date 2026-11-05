"""Conservative first-intersection reference for a single monotone flux-label branch."""
import numpy as np
from scipy.special import ndtr
from .broadening import _cdf_interval


def branch_power(edges, centers, widths, weights):
    """Deposit Gaussian displacement in the flux label, then map each bin to its first wall hit.

    Input labels and sigma share units. Edges are in wall traversal order. The caller
    must supply a verified one-to-one first-hit map. Both outside masses are returned
    separately, without a closure prescription at the separatrix or a second wall.
    """
    edges, centers, widths, weights = (np.asarray(v, dtype=float) for v in (edges, centers, widths, weights))
    if (edges.ndim!=1 or len(edges)<2 or centers.ndim!=1 or not len(centers)
            or widths.shape!=centers.shape or weights.shape!=centers.shape
            or not all(np.isfinite(x).all() for x in (edges,centers,widths,weights))
            or np.any(widths<=0) or np.any(weights<0) or not np.isfinite(weights.sum()) or weights.sum()<=0):
        raise ValueError('finite Gaussian donors, positive widths and positive total weight required')
    reverse = np.all(np.diff(edges)<0)
    if reverse:
        edges = edges[::-1]
    if np.any(np.diff(edges)<=0):
        raise ValueError('ambiguous first-hit branch: flux-coordinate edges must be strictly monotone')
    power=np.zeros(len(edges)-1)
    below=above=0.
    for x,sigma,w in zip(centers,widths,weights):
        z=(edges-x)/sigma
        power += w*_cdf_interval(z[:-1],z[1:])
        below += w*ndtr(z[0])
        above += w*ndtr(-z[-1])
    error=abs(power.sum()+below+above-weights.sum())/weights.sum()
    if error>2e-12 or np.any(power<0):
        raise FloatingPointError('flux-label ledger did not close')
    return dict(cell_power=power[::-1] if reverse else power,below=float(below),above=float(above),
                prescribed=float(weights.sum()),balance_relative=float(error))

