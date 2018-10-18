%% CORRER EN SERIE
es=1;
while es > 15.27366971818187e-20
x(1)=0;
while x(1) < x(4)
    scriptWT
end
res=errores(x);
es=sum(sum((data-res).^2));
end
prueba