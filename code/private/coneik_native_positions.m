function [Q,codes,native_precision,precision_unreachable,precision_dps]=coneik_native_positions(gamma,alpha)
% Keep the binary64 MEX ABI intact; attempt all unresolved rows in one call.
    position_backend=coneik_native_backend('position');
    [Q,codes]=position_backend(gamma,alpha);
    native_precision=false(size(codes));precision_unreachable=false(size(codes));
    precision_dps=nan(size(codes));
    rows=find(codes==-2);
    if isempty(rows),return;end
    precision_backend=coneik_native_backend('precision');
    if isempty(precision_backend),return;end
    [P,pc,pdps]=precision_backend(gamma(rows,:),alpha);
    solved=pc==0|pc==-1;
    Q(rows(solved),:)=P(solved,:);codes(rows(solved))=pc(solved);
    native_precision(rows(pc==0))=true;
    precision_unreachable(rows(pc==-1))=true;
    precision_dps(rows(solved))=pdps(solved);
end
