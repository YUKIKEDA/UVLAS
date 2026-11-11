"""Internal accepted-state reference; exact modal flow is a restart oracle only."""
import copy
import hashlib
import json
import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                     separators=(',', ':')).encode()).hexdigest()


class Reference:
    def __init__(self, mass, stiffness, modes, omega2, model_hash):
        self.M, self.K, self.V, self.w = mass, stiffness, modes, np.sqrt(omega2)
        self.identity = dict(model=model_hash, method='exact-modal-reference-v1',
                             dofs=len(mass), dtype='float64', history_version=1)

    def state(self, x, v, time=0., step=0, dt=0.):
        out = dict(self.identity, x=np.asarray(x).tolist(), v=np.asarray(v).tolist(),
                   time=float(time), step=int(step), cursor=int(step), dt=float(dt),
                   accepted=True, history={})
        self.validate(out)
        return out

    def validate(self, s):
        if any(s.get(k) != v for k, v in self.identity.items()):
            raise ValueError('identity')
        if s.get('accepted') is not True or s.get('history') != {}:
            raise ValueError('accepted/history')
        if type(s.get('step')) is not int or s['step'] < 0 or s.get('cursor') != s['step']:
            raise ValueError('generation/cursor')
        if not np.isfinite([s.get('time', np.nan), s.get('dt', np.nan)]).all() or min(s['time'], s['dt']) < 0:
            raise ValueError('time')
        for key in ['x', 'v']:
            a = np.asarray(s.get(key), dtype=float)
            if a.shape != (self.identity['dofs'],) or not np.isfinite(a).all():
                raise ValueError('array')

    def trial(self, accepted, dt):
        self.validate(accepted)
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('step size')
        a = self.V.T @ self.M @ np.asarray(accepted['x'])
        b = self.V.T @ self.M @ np.asarray(accepted['v'])
        c, s = np.cos(self.w*dt), np.sin(self.w*dt)
        x = self.V @ (c*a+s*b/self.w)
        v = self.V @ (-self.w*s*a+c*b)
        return dict(parent=digest(accepted), candidate=self.state(x, v, accepted['time']+dt,
                    accepted['step']+1, dt), used=False)

    def commit(self, accepted, trial):
        if trial['used'] or trial['parent'] != digest(accepted):
            raise ValueError('stale/double commit')
        self.validate(trial['candidate'])
        if trial['candidate']['step'] != accepted['step']+1:
            raise ValueError('generation')
        trial['used'] = True
        return copy.deepcopy(trial['candidate'])

    def save(self, state):
        self.validate(state)
        return json.dumps(dict(payload=state, sha256=digest(state)), allow_nan=False)

    def load(self, text):
        wrapper = json.loads(text)
        if digest(wrapper['payload']) != wrapper['sha256']:
            raise ValueError('checksum')
        self.validate(wrapper['payload'])
        return copy.deepcopy(wrapper['payload'])
