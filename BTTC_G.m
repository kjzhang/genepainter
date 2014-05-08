function g = BTTC_G(image, P, P1, P2, P3)
    %{
    c1 = get_c(image, P1);
    c2 = get_c(image, P2);
    c3 = get_c(image, P3);
    a = get_alpha(P, P1, P2, P3);
    b = get_beta(P, P1, P2, P3);
    g = c1 + a * (c2 - c1) + b * (c3 - c1);
    %}

    g = G_approx(P1, P2, P3, P, image);
end

%{
function a = get_alpha(P, P1, P2, P3)
    x = P(1);
    y = P(2);
    
    x1 = P1(1);
    y1 = P1(2);
    
    x2 = P2(1);
    y2 = P2(2);
    
    x3 = P3(1);
    y3 = P3(2);
    
    a = ((x - x1) * (y3 - y1) - (y - y1) * (x3 - y1)) / ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))
end

function b = get_beta(P, P1, P2, P3)
    x = P(1);
    y = P(2);
    
    x1 = P1(1);
    y1 = P1(2);
    
    x2 = P2(1);
    y2 = P2(2);
    
    x3 = P3(1);
    y3 = P3(2);
    
    b = ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) / ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))
end

function c = get_c(image, p)
    c = image(p(1), p(2));
end
%}

function a = alpha(P1,P2,P3,p)
    a = (p(1)-P1(1))*(P3(2)-P1(2)) - (p(2)-P1(2))*(P3(1)-P1(1));
    a = a/((P2(1)-P1(1))*(P3(2)-P1(2)) - (P2(2)-P1(2))*(P3(1)-P1(1)));
end

function b = beta(P1,P2,P3,p)
    b = (P2(1)-P1(1))*(p(2)-P1(2)) - (P2(2)-P1(2))*(p(1)-P1(1));
    b = b/((P2(1)-P1(1))*(P3(2)-P1(2)) - (P2(2)-P1(2))*(P3(1)-P1(1)));
end

function Gxy = G_approx(P1,P2,P3,p,F)
    c1 = F(P1(1),P1(2),:);
    c2 = F(P2(1),P2(2),:);
    c3 = F(P3(1),P3(2),:);
    a = alpha(P1,P2,P3,p);
    b = beta(P1,P2,P3,p);
    Gxy = c1 + a*(c2-c1) + b*(c3-c1);
end