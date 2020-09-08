%% linear unbiased MMSE detector (LMMSE)
function [idxhat,bithat] = LMMSE(par,H,y,N0)
  W = (H'*H+(N0/par.Es)*eye(par.U))\(H');
  xhat = W*y;
  G = real(diag(W*H));
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end