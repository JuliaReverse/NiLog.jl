using NiLogLikeNumbers
using Test

@testset "einsum" begin
	N = 3
	a = Tropical.(randn(N, N))
	b = Tropical.(randn(N, N))
	c = Tropical.(zeros(N, N))
	einsum!(((1,2), (2,3)), (a, b), (1,3), c)
	@test c ≈ a*b
end
