import copy
import unittest
import tempfile
from pathlib import Path
from qualification_selection import ranked_improvements, passes_improvement


def validation():
    return dict(complete=True, verdict=dict(resolution_screen_passed=True, improvement_confirmed=True),
        levels=[dict(method=m, tubes=n, rows={k:dict(audit_passed=True, feasible=True,
            peak_MW_m2=p) for k,p in [('baseline',2.7),('selected',1.8)]},
            geometry_audit=dict(passed_independent=True,geometry_screen=dict(passed=True)))
            for m,n in [('cubic',320),('cubic',640),('cubic',1280),('quintic',640)]],
        wall_only_refinement=dict(geometry_audit=dict(passed_independent=True),
            rows={k:dict(audit_passed=True,feasible=True) for k in ['baseline','selected']}))


class SelectionTests(unittest.TestCase):
    def test_fine_order_not_coarse_order(self):
        rows=[dict(identity='a',description=dict(feasible=True,geometry_screen=dict(passed=True),peak_MW_m2=2.0,peak_ratio=.8)),
              dict(identity='b',description=dict(feasible=True,geometry_screen=dict(passed=True),peak_MW_m2=1.8,peak_ratio=.7))]
        self.assertEqual([r['identity'] for r in ranked_improvements(rows)],['b','a'])
        rows[1]['description']['geometry_screen']['passed']=False
        self.assertEqual([r['identity'] for r in ranked_improvements(rows)],['a'])

    def test_numeric_pass_without_improvement_not_promoted(self):
        v=validation();self.assertTrue(passes_improvement(v))
        v['verdict']['improvement_confirmed']=False
        self.assertFalse(passes_improvement(v))

    def test_all_checks_required(self):
        for mutate in [lambda v:v.update(complete=False),
                       lambda v:v['levels'].pop(),
                       lambda v:v['levels'][2]['rows']['selected'].update(peak_MW_m2=3.0),
                       lambda v:v['levels'][2]['rows']['selected'].update(peak_MW_m2=float('nan')),
                       lambda v:v['levels'][0]['geometry_audit'].update(passed_independent=False),
                       lambda v:v['wall_only_refinement']['geometry_audit'].update(passed_independent=False),
                       lambda v:v['levels'][3]['rows']['baseline'].update(audit_passed=False)]:
            v=validation();mutate(v);self.assertFalse(passes_improvement(v))

    def test_failed_best_falls_back_to_next_passing(self):
        invalid=validation();invalid['verdict']['improvement_confirmed']=False
        attempts=[('best',invalid),('next',validation()),('worse',validation())]
        selected=next((identity for identity,v in attempts if passes_improvement(v)),None)
        self.assertEqual(selected,'next')

    def test_cached_failure_and_success_records_are_checksummed(self):
        from qualify_omp_stages import checked_read,checked_save,read,save
        with tempfile.TemporaryDirectory(dir='D:/DivertorOptimization') as temp:
            p=Path(temp)/'validation.json';payload=validation();checked_save(p,payload)
            self.assertEqual(checked_read(p),payload)
            damaged=read(p);damaged['payload']['verdict']['improvement_confirmed']=False;save(p,damaged)
            with self.assertRaisesRegex(ValueError,'checksum'):checked_read(p)


if __name__=='__main__': unittest.main()
