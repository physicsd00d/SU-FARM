
from scipy.stats import f

# # Reference values, should give Pmid = 0.759078400489
# Nf = 20   # Num Failures
# Nt = 26   # Num Total launches
# conf = 0.6

Nf = 1
Nt = 1e6

conf = 0.6
alpha = (1-conf)/2

vu1 = 2*Nf + 2
vu2 = 2*(Nt - Nf)

vd1 = 2*(1+Nt-Nf)
vd2 = 2*Nf

fup = f.isf(alpha,vu1,vu2)
flow = f.isf(alpha,vd1,vd2)

Pup = (Nf+1.)*fup/(Nt - Nf + (Nf+1.)*fup)
Plow = Nf/(Nf + (Nt - Nf +1)*flow)
Pmid = (Pup+Plow)/2.

print "Pmid = {0}".format(Pmid)
