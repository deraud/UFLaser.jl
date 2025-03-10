const a = 2
const b = 5

Base.@kwdef mutable struct test
    c = 5
end

function adding!(t::test)
    t.c = a + t.c +1
end
testt = test()
adding!(testt)
adding!(testt)
println(testt)