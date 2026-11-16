# REPL displays

_sty(io, s, color; bold=false) = get(io, :color, false) ?
    string(bold ? "\e[1m" : "", "\e[38;5;$(color)m", s, "\e[0m") : string(s)
_pfmt(p) = isnan(p) ? "" : p < 0.001 ? "<0.001" : @sprintf("%.3f", p)
_stars(p) = isnan(p) ? "" : p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : p < 0.1 ? "." : ""
_f(x, fmt) = isnan(x) ? "" : @sprintf(fmt, x)

Base.show(io::IO, r::PermanovaResult) = print(io, "PermanovaResult(", join(r.terms[1:end-2], ", "), ")")

const _BYTEXT = Dict(:terms => "terms added sequentially (first to last)",
                     :margin => "marginal effects of terms", :model => "whole model")

function Base.show(io::IO, ::MIME"text/plain", r::PermanovaResult)
    println(io, _sty(io, "PERMANOVA", 117; bold=true), _sty(io, " · ", 240),
            "$(r.n) observations, $(r.permutations) permutations, $(_BYTEXT[r.by])")
    println(io, _sty(io, "  model: ~ " * r.formula, 245))
    w = max(8, maximum(length, r.terms))
    println(io, _sty(io, "  " * rpad("", w) * "    Df   SumOfSqs       R²        F   Pr(>F)", 245))
    for i in eachindex(r.terms)
        st = _stars(r.p[i])
        line = "  " * rpad(r.terms[i], w) * lpad(r.df[i], 6) * lpad(_f(r.sumsq[i], "%.4f"), 11) *
               lpad(_f(r.r2[i], "%.4f"), 9) * lpad(_f(r.F[i], "%.3f"), 9) * lpad(_pfmt(r.p[i]), 9)
        print(io, line, " ", _sty(io, st, 221))
        i < length(r.terms) && println(io)
    end
    any(<(0.1), filter(!isnan, r.p)) &&
        print(io, "\n", _sty(io, "  signif.: *** <0.001  ** <0.01  * <0.05  . <0.1", 242))
end

Base.show(io::IO, r::PairwisePermanova) = print(io, "PairwisePermanova($(r.variable), $(length(r.p)) pairs)")

function Base.show(io::IO, ::MIME"text/plain", r::PairwisePermanova)
    println(io, _sty(io, "Pairwise PERMANOVA", 117; bold=true), _sty(io, " · ", 240),
            "$(r.variable), $(length(r.p)) pairs, p adjusted: $(r.adjust)")
    names = [a * " – " * b for (a, b) in zip(r.group1, r.group2)]
    w = max(4, maximum(length, names))
    print(io, _sty(io, "  " * rpad("pair", w) * "     n        F       R²        p    p.adj", 245))
    for i in eachindex(names)
        st = _stars(r.p_adjusted[i])
        print(io, "\n  ", rpad(names[i], w), lpad(r.n[i], 6), lpad(_f(r.F[i], "%.3f"), 9),
              lpad(_f(r.r2[i], "%.4f"), 9), lpad(_pfmt(r.p[i]), 9), lpad(_pfmt(r.p_adjusted[i]), 9), " ",
              _sty(io, st, 221))
    end
end
