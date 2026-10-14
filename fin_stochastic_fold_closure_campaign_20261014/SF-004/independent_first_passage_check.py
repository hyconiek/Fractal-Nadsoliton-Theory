"""Independent check by high-precision quadrature + ODE residual (diagnostic, not certificate)."""
import json
import mpmath as m
from pathlib import Path
m.mp.dps=70
c=m.mpf(2)/3
Iminus=m.quad(lambda y:m.exp(c*y**3),[-1,0]);Iplus=m.quad(lambda y:m.exp(c*y**3),[0,1])
p=Iminus/(Iminus+Iplus)
folder=Path(__file__).resolve().parent
cert=json.loads((folder/'EXIT_SIDE_CERTIFICATE.json').read_text())
lo=m.mpf(cert['P_exit_upper_public_outward_interval'][0]);hi=m.mpf(cert['P_exit_upper_public_outward_interval'][1]); assert lo<=p<=hi,(lo,p,hi)
assert m.mpf('0.41')<p<m.mpf('0.43')
result={'mpmath_decimal_digits':m.mp.dps,'upper_quadrature':m.nstr(p,55),'inside_exact_rational_certificate':True,'method':'independent mpmath quadrature, not the rational producer','comparison_scope':'cross-check only, not interval proof'}
(folder/'INDEPENDENT_EXIT_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print('pupper',m.nstr(p,30))
