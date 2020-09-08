
%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
  xhat = H\y;    
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end