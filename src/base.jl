export igemm!, igemv!, isum
export maxloc

const LogLikeNumber{T} = Union{ULogarithmic{T}, Tropical{T}}

const TropicalG{T,TG} = Tropical{GVar{T,TG}}

@i @inline function :(*=)(+)(z::T, x::T, y::T) where T<:Tropical
    @invcheckoff if (content(x) > content(y), ~)
        content(z) += identity(content(x))
    else
        content(z) += identity(content(y))
    end
end

@i @inline function (:*=(identity))(x::T, y::T) where T<:LogLikeNumber
    content(x) += identity(content(y))
end

@i @inline function (:*=(*))(out!::T, x::T, y::T) where T<:LogLikeNumber
    content(out!) += content(x) + content(y)
end

export tropical_muladd
# branch should be initialized to false.
@i @inline function tropical_muladd(out!::T, x::T, y::T, branch::Bool) where T<:Tropical
	x *= identity(y)
	@invcheckoff if (out! < x, branch)
		FLIP(branch)
		NiLang.SWAP(out!, x)
	end
end

maxloc(v::AbstractVector) = findmax(v)[2]

# Happy Pirate!!!!!
ULogarithmic{T}(gv::T) where T<:Real = exp(ULogarithmic{T}, gv)
ULogarithmic(gv::T) where T<:Real = exp(ULogarithmic, gv)

import TropicalNumbers: content
@fieldview content(x::ULogarithmic) = x.log
NiLang.chfield(x::Tropical{T}, ::typeof(content), val::T) where T = Tropical{T}(val)

# AD wrappers
for T in [:Tropical, :ULogarithmic]
    @eval NiLang.AD.GVar(x::$T) = $T(GVar(content(x), zero(content(x))))
    @eval (_::Type{Inv{$T}})(x::$T) = content(x)
    @eval NiLang.AD.grad(x::$T{<:GVar}) = $T(grad(content(x)))
    @eval (_::Type{Inv{GVar}})(x::$T{<:GVar}) = $T((~GVar)(content(x)))

    @eval Base.one(x::$T{GVar{T,GT}}) where {T, GT} = one($T{GVar{T,GT}})
    @eval Base.one(::Type{$T{GVar{T,GT}}}) where {T,GT} = $T(GVar(zero(T), zero(GT)))
    @eval Base.zero(x::$T{GVar{T,GT}}) where {T,GT} =zero($T{GVar{T,GT}})
    @eval Base.zero(::Type{$T{GVar{T,T}}}) where T = GVar(zero(T))
end

function NiLang.loaddata(::Type{Array{LogLikeNumber{GVar{T,T}},N}}, data::Array{LogLikeNumber{T},N}) where {T,N}
    GVar.(data)
end
import NiLang.NiLangCore: deanc

function deanc(x::T, v::T) where T<:LogLikeNumber
    x === v || deanc(content(x), content(v))
end
