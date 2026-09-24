# The certificate is algebraic: Delta U=0 and I>0 imply c=0.
gamma=3.0; I=2.5; deltaU=0.0
c=deltaU/(gamma*I)
assert c==0.0
for deltaU in [1,-1]:
 assert deltaU/(gamma*I)!=0
print('PASS DYN-004 no-moving-branch certificate')
