# ============================
# 对称 3×3 Jacobi 特征分解
# ============================
# 输入：上三角 6 个分量 + sweep 次数
# 输出：按降序排列的特征值 (b1,b2,b3) 以及特征向量矩阵 Q 的 9 个标量（按列存储）
@inline Base.@propagate_inbounds function sym_eig3x3_jacobi(
    A11::T2, A22::T2, A33::T2, A12::T2, A13::T2, A23::T2, sweeps::Integer=6
) where T2
    a11 = T2(A11); a22 = T2(A22); a33 = T2(A33)
    a12 = T2(A12); a13 = T2(A13); a23 = T2(A23)

    # Q 初始化为单位阵
    q11 = one(T2) ; q12 = zero(T2); q13 = zero(T2)
    q21 = zero(T2); q22 = one(T2) ; q23 = zero(T2)
    q31 = zero(T2); q32 = zero(T2); q33 = one(T2)

    @inbounds for _ = 1:sweeps
        # 轮 1：消 (1,2)
        a11,a22,a12, a13,a23, q11,q21,q31,q12,q22,q32 =
            jacobi_pq!(a11,a22,a12, a13,a23, q11,q21,q31,q12,q22,q32)
        # 轮 2：消 (1,3)
        a11,a33,a13, a12,a23, q11,q21,q31,q13,q23,q33 =
            jacobi_pq!(a11,a33,a13, a12,a23, q11,q21,q31,q13,q23,q33)
        # 轮 3：消 (2,3)
        a22,a33,a23, a12,a13, q12,q22,q32,q13,q23,q33 =
            jacobi_pq!(a22,a33,a23, a12,a13, q12,q22,q32,q13,q23,q33)
    end

    # 现在近似对角：a11,a22,a33；Q 的列为特征向量
    # 按降序排序，并同步交换 Q 的列
    b1 = a11; b2 = a22; b3 = a33
    if b1 < b2
        b1,b2 = b2,b1
        q11,q12 = q12,q11; q21,q22 = q22,q21; q31,q32 = q32,q31
    end
    if b1 < b3
        b1,b3 = b3,b1
        q11,q13 = q13,q11; q21,q23 = q23,q21; q31,q33 = q33,q31
    end
    if b2 < b3
        b2,b3 = b3,b2
        q12,q13 = q13,q12; q22,q23 = q23,q22; q32,q33 = q33,q32
    end

    return b1,b2,b3, q11,q12,q13, q21,q22,q23, q31,q32,q33
end

# 单次 Jacobi 旋转：消去 (p,q) 元（这里参数已按位展开）
@inline Base.@propagate_inbounds function jacobi_pq!(
    app::T2, aqq::T2, apq::T2, apm::T2, aqm::T2,
    q1p::T2, q2p::T2, q3p::T2, q1q::T2, q2q::T2, q3q::T2
) where T2
    thresh = eps(T2) * max(abs(app) + abs(aqq), one(T2))
    if abs(apq) > thresh
        τ = (aqq - app) / (T2(2)*apq)
        t = sign(τ) / (abs(τ) + sqrt(T2(1) + τ*τ))
        c = 1/(sqrt(T2(1) + t*t))
        s = t*c

        app2 = app - t*apq
        aqq2 = aqq + t*apq
        apq2 = T2(0)

        apm2 = c*apm - s*aqm
        aqm2 = s*apm + c*aqm

        q1p2 = c*q1p - s*q1q;  q1q2 = s*q1p + c*q1q
        q2p2 = c*q2p - s*q2q;  q2q2 = s*q2p + c*q2q
        q3p2 = c*q3p - s*q3q;  q3q2 = s*q3p + c*q3q

        return app2, aqq2, apq2, apm2, aqm2, q1p2, q2p2, q3p2, q1q2, q2q2, q3q2
    else
        return app, aqq, apq, apm, aqm, q1p, q2p, q3p, q1q, q2q, q3q
    end
end