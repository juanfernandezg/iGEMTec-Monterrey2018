%% CORRER EN SERIE
es=1;
while es > 1.5e-19
x(1)=0;
while x(1) < x(4)
    scriptRT
end
res=errores(x);
es=sum(sum((data-res).^2));
end
prueba