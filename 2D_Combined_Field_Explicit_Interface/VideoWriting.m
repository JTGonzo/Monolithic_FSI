aviobj = avifile('example.avi','compression','None');
% nFrames = 20;
% 
% % Preallocate movie structure.
% mov(1:nFrames) = struct('cdata', [],...
%                         'colormap', []);
% 
% % Create movie.
% Z = peaks; surf(Z); 
% axis tight
% set(gca,'nextplot','replacechildren');
for k = 1:nFrames 
   surf(sin(2*pi*k/20)*Z,Z)
   aviobj = addframe(aviobj,gcf);
   drawnow 
%    mov(k) = getframe(gcf);
end

% Create AVI file.
viobj = close(aviobj);

