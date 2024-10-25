# Julia code to calculate prolongation stencil coefficients

using LinearAlgebra

# General assumptions:
# - We label vertices. Cells (cell centred) are offset by 1/2 to the right.
# - We assume the refinement factor is 2.
# - Fine points have an offset of either 0 or 1 with respect to the coarse grid.
# - We only calculate vertex centred stencils for `offset = 1`. The
#   case `offset = 0` is trivial.
# - We only calculate cell centred stencils for `offset = 1`. The case
#   `offset = 0` follows by symmetry, and is handled when the stencil
#   is evaluated.

# Vertex-centred layout:
#     coarse position:    -2      -1       0      +1      +2
#                          X       X       X       X       X
#                          x   x   x   x   x   x   x   x   x
#     fine offset:                         0   1

# Cell-centred layout:
#     coarse position:        -2      -1       0      +1      +2
#                              O       O       O       O       O
#                                o   o   o   o   o   o   o   o
#     fine offset:                           0   1

# For polynomial interpolation we do the following:
# - Set up a generic polynomial of the respective order `p`, i.e. with `p+1` coefficients
# - Choose `p+1` coarse grid points, and then choose the polynomial
#   coeefficiens so that it agrees with the coarse grid values at
#   these coarse grid points
# - Evaluate the polynomial at the desired point

const BigRat = Rational{BigInt}

function stencil_vc_poly(p::Integer, x)
    # Obviously
    @assert p >= 0

    # For symmetry we require an odd `p`
    @assert isodd(p)

    imin = -div(p - 1, 2)
    imax = +div(p + 1, 2)
    @assert imax - imin == p

    return stencil_vc_poly(imin, imax, x)
end

function stencil_vc_poly(imin::Integer, imax::Integer, x0::Real)
    # Obviously
    p = imax - imin
    @assert p >= 0

    x = BigRat(x0)

    # Prevent extrapolation
    imin <= x <= imax || @show imin x imax
    @assert imin <= x <= imax

    # Coarse grid points
    # We choose integer points for convenience. We could scale by a
    # grid spacing `h` but that would just cancel out.
    xs = imin:imax
    @assert length(xs) == p + 1

    # We want interpolation to be exact for functions of orders from
    # `0` to `p`.
    # This is the order of the function that we represent via a `p`-th
    # order polynomial, so its order can be anything from `0` to `p`.

    # Our test functions (polynomials)
    testfun(x, q::Integer) = x^q

    # Determine polynomial coefficients `cs`
    # - We are solving a linear system `A cs = b`
    # - Each row of the system corresponds to one condition.
    # - Each conditions states that, when we have a function of order
    #   `q`, the interpolation must be exact.
    A = fill(big(1 // 0), p + 1, p + 1)
    b = fill(big(1 // 0), p + 1)
    for q in 0:p
        # Function values at the grid points
        ys = testfun.(xs, q)
        # Evaluate at `x`
        y = testfun(x, q)
        # Condition: With the inputs `ys`, the result must be `y`
        for i in 0:p
            A[q + 1, i + 1] = ys[i + 1]
        end
        b[q + 1] = y
    end
    # Ensure all elements have been set
    @assert all(A .!= 1 // 0)
    @assert all(b .!= 1 // 0)

    cs = A \ b

    # Test the conditions
    for q in 0:p
        # Function values at the grid points
        ys = testfun.(xs, q)
        # Evaluate at `x`
        y = testfun(x, q)
        # Check evaluated point
        @assert y == sum(cs .* ys)
    end

    return cs::AbstractVector{BigRat}
end

function stencil_cc_poly(p::Integer, x::Real)
    # Obviously
    @assert p >= 0

    # Choosing `imin` and `imax`:
    if x < 0
        # Lean left
        #     p    xs
        #     0    -1:-1
        #     1    -1:+0
        #     2    -2:+0
        #     3    -2:+1
        imin = -fld(p, 2) - 1
        imax = +cld(p, 2) - 1
    elseif x > 0
        # Lean right
        #     p    xs
        #     0    -0:+0
        #     1    -0:+1
        #     2    -1:+1
        #     3    -1:+2
        imin = -fld(p, 2)
        imax = +cld(p, 2)
    else
        @assert iseven(p)
        imin = -div(p, 2)
        imax = +div(p, 2)
    end
    @assert imax - imin == p

    return stencil_cc_poly(imin, imax, x)
end

function stencil_cc_poly(imin::Integer, imax::Integer, x0::Real)
    # Obviously
    p = imax - imin
    @assert p >= 0

    x = BigRat(x0)

    # Prevent extrapolation
    @assert imin <= x <= imax + 1

    # Coarse grid points
    # We choose integer points for convenience. We could scale by a
    # grid spacing `h` but that would just cancel out.
    xs = imin:imax
    @assert length(xs) == p + 1

    # We want interpolation to be exact for functions of orders from
    # `0` to `p`.
    # This is the order of the function that we represent via a `p`-th
    # order polynomial, so its order can be anything from `0` to `p`.

    # Our test functions (polynomials)
    testfun(x, q::Integer) = x^q

    # Determine polynomial coefficients `cs`
    # - We are solving a linear system `A cs = b`
    # - Each row of the system corresponds to one condition.
    # - Each conditions states that, when we have a function of order
    #   `q`, the interpolation must be exact.
    A = fill(big(1 // 0), p + 1, p + 1)
    b = fill(big(1 // 0), p + 1)
    for q in 0:p
        # Function values at the grid points
        ys = testfun.(xs, q)
        # Evaluate at `x`
        y = testfun(x, q)
        # Condition: With the inputs `ys`, the result must be `y`
        for i in 0:p
            A[q + 1, i + 1] = ys[i + 1]
        end
        b[q + 1] = y
    end
    # Ensure all elements have been set
    @assert all(A .!= 1 // 0)
    @assert all(b .!= 1 // 0)

    cs = A \ b

    # Test the conditions
    for q in 0:p
        # Function values at the grid points
        ys = testfun.(xs, q)
        # Evaluate at `x`
        y = testfun(x, q)
        # Check evaluated point
        @assert y == sum(cs .* ys)
    end

    return cs::AbstractVector{BigRat}
end

function stencil_cc_cons(p::Integer, x::Real)
    # Obviously
    @assert p >= 0

    # Choosing `imin` and `imax`:
    if x < 0
        # Lean left
        #     p    xs
        #     0    -0:+0
        #     1    -1:+0
        #     2    -1:+1
        #     3    -2:+2
        imin = -cld(p, 2)
        imax = +fld(p, 2)
    elseif x > 0
        # Lean right
        #     p    xs
        #     0    -0:+0
        #     1    -0:+1
        #     2    -1:+1
        #     3    -1:+2
        imin = -fld(p, 2)
        imax = +cld(p, 2)
    else
        @assert iseven(p)
        imin = -div(p, 2)
        imax = +div(p, 2)
    end

    return stencil_cc_cons(imin, imax, x)
end

function stencil_cc_cons(imin::Integer, imax::Integer, x0::Real)
    # Obviously
    p = imax - imin
    @assert p >= 0

    x = BigRat(x0)

    # Prevent extrapolation
    @assert imin <= x <= imax + 1

    # We implement conservative cell-centred interpolation by falling
    # back to polynomial vertex-centred interpolation.
    # - First calculate the antiderivative, switching to vertex
    #   centring, and increasing the order by one, and the stencil
    #   size by one to the right
    # - The interpolate there polynomially
    # - Then calculate the derivative, switching back to cell centring

    # The antiderivative operator
    antideriv = [col < row ? 1 : 0 for row in 0:(p + 1), col in 0:p]

    # The derivative operator
    deriv = [col == row ? -1 : col == row + 1 ? +1 : 0 for row in 0:p, col in 0:(p + 1)]

    @assert deriv * antideriv == I

    # This is the derivative operator, written out explicitly
    x_left = x - 1 // 4
    x_right = x + 1 // 4
    cs_poly_left = stencil_vc_poly(imin, imax + 1, x_left)
    cs_poly_right = stencil_vc_poly(imin, imax + 1, x_right)

    # Multiply by the refinement factor `2`, i.e. divide by the ratio
    # of the grid spacings
    cs_poly_matrix = reshape((cs_poly_right - cs_poly_left) / (x_right - x_left), 1, :)
    cs_matrix = cs_poly_matrix * antideriv
    cs = reshape(cs_matrix, :)

    return cs::AbstractVector{BigRat}
end

function main()
    display([(p=p, cs=stencil_vc_poly(p, 1 // 2)) for p in 1:2:5])
    display([(p=p, cs=stencil_cc_poly(p, -1 // 4)) for p in 0:5])
    display([(p=p, cs=stencil_cc_poly(p, +1 // 4)) for p in 0:5])
    display([(p=p, s=s, cs=stencil_cc_cons(-p + s - 1, -1 + s, -1 // 4)) for p in 0:5 for s in 0:p])
    display([(p=p, s=s, cs=stencil_cc_cons(-p + s, 0 + s, +1 // 4)) for p in 0:5 for s in 0:p])
    nothing
end

nothing
