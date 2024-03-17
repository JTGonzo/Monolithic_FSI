a1=round(19.5/Sol.delt);
a2=round(20/Sol.delt);
figure;
subplot(2,1,1); plot((a1:a2)*Sol.delt,Sol.ad(a1:a2,1))
title('x')
grid on;

subplot(2,1,2); plot((a1:a2)*Sol.delt,Sol.ad(a1:a2,2))
title('y')
grid on;

figure;
subplot(2,1,1); plot((a1:a2)*Sol.delt,-Sol.ldc(a1:a2,2)*1e3)
title('lift')
grid on;

subplot(2,1,2); plot((a1:a2)*Sol.delt,-Sol.ldc(a1:a2,1)*1e3)
title('drag')
grid on;

