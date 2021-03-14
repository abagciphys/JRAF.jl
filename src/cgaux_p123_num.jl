###########################################################################################################
###########################################################################################################
############################ OVERLAP LIKE G VİA CUBA NUMERICAL INTEGRATION ################################
###########################################################################################################
###########################################################################################################
export CuhreAuxiliaryG, VegasAuxiliaryG, SuaveAuxiliaryG
###########################################################################################################
function CuhreAuxiliaryG(n1, q, n2, n3, p1, p2, p3, acc)
    coef = ((p1^(n1))//Gamma(n1+1))
    function integrand(x, f)

        f[1] = ( (((1/(1 - x[1]))*(2*x[2] - 1))^q)*
                 (((1/(1 - x[1])) + (2*x[2] - 1))^n2)*
                 (((1/(1 - x[1])) - (2*x[2] - 1))^n3)*
                 exp(-p2*(1/(1 - x[1])) - p3*(2*x[2] - 1))*
                 (2/((1 - x[1])^2))
        )

    end
    res, err = cuhre(integrand, 2, 1, atol=acc, rtol=acc)
    result = NO(coef) * res

    return result
end

function VegasAuxiliaryG(n1, q, n2, n3, p1, p2, p3, acc)
    coef = ((p1^(n1))//Gamma(n1+1))
    function integrand(x, f)

        f[1] = ( (((1/(1 - x[1]))*(2*x[2] - 1))^q)*
                 (((1/(1 - x[1])) + (2*x[2] - 1))^n2)*
                 (((1/(1 - x[1])) - (2*x[2] - 1))^n3)*
                 exp(-p2*(1/(1 - x[1])) - p3*(2*x[2] - 1))*
                 (2/((1 - x[1])^2))
        )

    end
    res, err = vegas(integrand, 2, 1, atol=acc, rtol=acc)
    result = NO(coef) * res

    return result
end

function SuaveAuxiliaryG(n1, q, n2, n3, p1, p2, p3, acc)
    coef = ((p1^(n1))//Gamma(n1+1))
    function integrand(x, f)

        f[1] = ( (((1/(1 - x[1]))*(2*x[2] - 1))^q)*
                 (((1/(1 - x[1])) + (2*x[2] - 1))^n2)*
                 (((1/(1 - x[1])) - (2*x[2] - 1))^n3)*
                 exp(-p2*(1/(1 - x[1])) - p3*(2*x[2] - 1))*
                 (2/((1 - x[1])^2))
        )

    end
    res, err = suave(integrand, 2, 1, atol=acc, rtol=acc)
    result = NO(coef) * res

    return result
end


