% cleanLFP2(matfile)
matfile = 'LFP_P1';
load(matfile);
% clear the edges

epsilon = 1e-1;
[nx, ny] = size(LFP);

% FP2;
% for j = 31:50
%     for i = 1:28
%         if i<=1.54*(j-31)
%             LFP(i,j) = 1e-10;
%             FP(i,j) = 1e-10;
%             exx(i,j) = 1;
%             exy(i,j) = 1;
%             eyy(i,j) = 1;
%             theta(i,j) = 0;
%             
%         end
%     end
% end

%% FP1
% for j = 1:16
%     for i = 85:110
%         if i>=1.5625*j+86
%             LFP(i,j) = 1e-10;
%             FP(i,j) = 1e-10;
%             exx(i,j) = 1;
%             exy(i,j) = 1;
%             eyy(i,j) = 1;
%             theta(i,j) = 0;
%             
%         end
%     end
% end


figure(3);
th = LFP+FP;
surf(th);view([0,0,1]);colormap;colorbar();
% save(matfile,'exx','eyy','exy','FP','HAADF','LFP','mask','theta','VADF','VBF');