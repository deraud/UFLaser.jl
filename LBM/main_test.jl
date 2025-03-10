Base.@kwdef mutable struct Mystruct
    a = 1
    b = 10
    c = 100
    d = 0
end

function additions!(m::Mystruct)
    m.d = m.a + m.b + m.d
    for i in 1:2
        m.d = m.d + 1
    end
end

function final(my::Mystruct)
    additions!(my)

end

mystruct = Mystruct()


final(mystruct)


println(mystruct.d)