function f = F_score(w,wt)
w = max(w,0);
beta = 2;
tp = sign(w)'*sign(wt);
fn = not(sign(w))'*sign(wt);
fp = sign(w)'*not(sign(wt));

f = ((1+beta^2)*tp)/((1+beta^2)*tp + (beta^2)*fn + fp);

end