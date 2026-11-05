"""Replay the published 28-point grid's interior-only selection from saved scores."""
import argparse
import json
from pathlib import Path


def interior_result(result):
    weights = [(6-i-j, i, j) for i in range(7) for j in range(7-i)]
    interior = {i for i, triple in enumerate(weights) if min(triple) > 0}

    def ranked(rows):
        if sorted(row['i'] for row in rows) != list(range(len(weights))):
            raise ValueError('expected the 28 original grid indices exactly once')
        return sorted((row for row in rows if row['eligible'] and row['i'] in interior),
                      key=lambda row: (row['peak'], row['i']))

    shortlist = {row['i'] for row in ranked(result['coarse'])[:3]}
    fine = ranked(result['fine'])
    chosen = next(row for row in fine if row['i'] in shortlist)
    return {
        'interior_posthoc_n': len(fine),
        'interior_chosen': chosen['i'],
        'interior_barycentric': list(weights[chosen['i']]),
        'interior_peak': chosen['peak'],
        'interior_finite_best': fine[0]['i'],
        'interior_excess': chosen['peak'] / fine[0]['peak'] - 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid', required=True, type=Path)
    parser.add_argument('--publication', required=True, type=Path,
                        help='directory containing corrected verified.json and REPORT.ko.md')
    args = parser.parse_args()
    summary = json.loads((args.publication / 'verified.json').read_text(encoding='utf-8'))
    report = (args.publication / 'REPORT.ko.md').read_text(encoding='utf-8')
    for row in summary['rows']:
        label = row['label']
        raw = json.loads((args.grid / label / 'result.json').read_text(encoding='utf-8'))
        actual = interior_result(raw)
        for field, value in actual.items():
            if row.get(field) != value:
                raise ValueError(f'{label} {field}: saved {row.get(field)!r}, recomputed {value!r}')
        table_row = (f"| {label} | {actual['interior_barycentric']} | {actual['interior_posthoc_n']} "
                     f"| {actual['interior_peak']:.4f} | {100*actual['interior_excess']:+.2f}% |")
        if table_row not in report.splitlines():
            raise ValueError(f'interior selection missing from report: {table_row}')
        print(table_row)


if __name__ == '__main__':
    main()
