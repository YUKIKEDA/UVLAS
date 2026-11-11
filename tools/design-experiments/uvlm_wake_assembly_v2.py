"""Frozen 2.2 numerical kernels, with stricter restart-policy validation."""
from uvlm_wake_assembly import *
from uvlm_wake_assembly import Model as BaseModel

class Model(BaseModel):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        object.__setattr__(self,'identity',digest('restart-policy-validated-v2',self.identity))
    def validate_state(self,state):
        result=super().validate_state(state)
        if state.initialization=='continued' and (state.epoch==0 or state.gamma_previous is None or state.last_h is None):raise ValueError('incomplete_accepted_circulation_history')
        for w in state.wakes:
            if state.max_intervals is not None and len(w.vertices)>state.max_intervals+1:raise ValueError('wake_exceeds_saved_window')
            for r in range(1,len(w.vertices)):
                row=w.row_ids[r]
                if not isinstance(row,tuple) or not row or row[0]!=w.sid or any(ident!=(row,j) for j,ident in enumerate(w.vertex_ids[r])):raise ValueError('invalid_material_row_ownership')
        return result
