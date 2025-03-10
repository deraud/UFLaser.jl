const a = 0
for i in 1:9
    a += lattice.f[1,1,i]
    print(a)
end