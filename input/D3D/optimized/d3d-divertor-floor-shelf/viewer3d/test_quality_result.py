import copy,hashlib,json,unittest
from quality_result import complete_report

def fixture():
    wall=[[1.,0.],[2.,0.]];digest=hashlib.sha256(json.dumps(wall,separators=(',',':')).encode()).hexdigest()
    summary={str(mm):{'all_wall_peak_q_over_P':1.,'peak_R_m':1.,'peak_Z_m':0.} for mm in [1,2,3,5,10]}
    patches=[{'component':'shelf','status':'direct_hit','angle_deg':1.6,'q':{str(mm):1. for mm in [1,2,3,5,10]}},
             {'component':'floor','status':'blocked_first_hit','angle_deg':None,'q':{str(mm):None for mm in [1,2,3,5,10]}}]
    models={k:{'shape_sha256':digest,'patches':copy.deepcopy(patches),'summary':copy.deepcopy(summary)} for k in ['baseline','candidate']}
    brief={k:{'peaks':[1.]*5,'locations':[[1.,0.]]*5,'unmapped_power':[0.]*5,'source_accounting_error':[0.]*5,
              'shelf_min_angle_deg':1.6,'component_power':{'floor':[1.]*5,'shelf':[1.]*5}} for k in models}
    surfaces={k:{'quadrature_levels':[{'integrals':{'floor':[1.]*5,'shelf':[1.]*5}}]*2} for k in models}
    report={'case_id':'case','shot_sha256':'shot','models':models,'levels':[{'models':brief}],
            'surfaces':surfaces,'forward_converged':True,'surface_integral_passed':True,'unresolved_target_length_m':0.}
    bundle={'id':'case','shot_sha256':'shot','models':{k:{'wall':wall} for k in models}}
    return report,bundle

class ResultTests(unittest.TestCase):
    def test_reference_angle_violation_does_not_invalidate_verified_candidate(self):
        r,b=fixture();r['models']['baseline']['patches'][0]['angle_deg']=1.3
        out=complete_report(r,b);self.assertTrue(out['verified']);self.assertFalse(out['baseline_constraint_feasible'])
    def test_compatible_report_passes(self):
        r,b=fixture();self.assertTrue(complete_report(r,b)['verified'])
    def test_stale_geometry_is_rejected(self):
        r,b=fixture();b['models']['candidate']['wall']=[[1.,0.],[2.,.001]]
        with self.assertRaisesRegex(ValueError,'Shape changed'):complete_report(r,b)
    def test_missing_connection_cannot_be_painted_zero(self):
        r,b=fixture();r['models']['candidate']['patches'][1]['q']['3']=0.
        with self.assertRaisesRegex(ValueError,'null heat'):complete_report(r,b)
    def test_shelf_violation_prevents_verification(self):
        r,b=fixture();r['models']['candidate']['patches'][0]['angle_deg']=1.49
        self.assertFalse(complete_report(r,b)['verified'])
    def test_independent_peak_disagreement_prevents_verification(self):
        r,b=fixture();r['models']['candidate']['summary']['3']['all_wall_peak_q_over_P']=1.02
        self.assertFalse(complete_report(r,b)['verified'])

if __name__=='__main__':unittest.main()
