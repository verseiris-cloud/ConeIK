function pp=coneik_boundary_spline(t,q,v0,vf)
%CONEIK_BOUNDARY_SPLINE Internal C2 cubic with one or two endpoint velocities.
% Empty velocity means not-a-knot on that side (n>=3). With two points,
% one velocity uses M1=M2, giving a quadratic; both give cubic Hermite.
% Physical seconds and nominal degrees. Numerical difficulty is unresolved.
    t=double(t(:));q=double(q(:));n=numel(t);
    if n<2||numel(q)~=n||any(~isfinite([t;q]))||any(diff(t)<=0)
        unresolved('Invalid internal interpolation data.');
    end
    if isempty(v0)&&isempty(vf),pp=spline(t,q);return;end
    h=diff(t);normal(h,true(size(h)),'interval');
    dq=diff(q);normal(dq,dq~=0,'position difference');
    delta=dq./h;normal(delta,dq~=0,'secant');
    A=spalloc(n,n,3*n);rhs=zeros(n,1);
    for i=2:n-1
        scale=max(h(i-1:i));ratio=h(i-1:i)/scale;
        normal(ratio,true(2,1),'interval ratio');
        A(i,i-1:i+1)=[ratio(1),2*sum(ratio),ratio(2)];
        difference=delta(i)-delta(i-1);normal(difference,difference~=0,'secant difference');
        rhs(i)=6*(difference/scale);normal(rhs(i),difference~=0,'interior equation');
    end
    if ~isempty(v0)
        A(1,1:2)=[2 1];difference=delta(1)-double(v0);
        rhs(1)=6*(difference/h(1));normal(rhs(1),difference~=0,'initial equation');
    elseif n==2
        A(1,:)=[1 -1];
    else
        scale=max(h(1:2));r=h(1:2)/scale;normal(r,true(2,1),'left knot ratio');
        A(1,1:3)=[-r(2),sum(r),-r(1)];
    end
    if ~isempty(vf)
        A(n,n-1:n)=[1 2];difference=double(vf)-delta(end);
        rhs(n)=6*(difference/h(end));normal(rhs(n),difference~=0,'final equation');
    elseif n==2
        A(n,:)=[1 -1];
    else
        scale=max(h(end-1:end));r=h(end-1:end)/scale;normal(r,true(2,1),'right knot ratio');
        A(n,n-2:n)=[r(2),-sum(r),r(1)];
    end
    if any(~isfinite(nonzeros(A)))||any(~isfinite(rhs)),unresolved('Nonfinite boundary system.');end
    try,M=A\rhs;catch,unresolved('Boundary linear solve failed.');end
    normal(M,M~=0,'second derivative solve');
    residual=A*M-rhs;
    scale=abs(A)*abs(M)+abs(rhs);
    if any(~isfinite([residual;scale]))||any(abs(residual)>1e-11*max(1,scale))
        unresolved('Boundary linear solve residual is unreliable.');
    end
    change=diff(M);a=(change./h)/6;normal(a,change~=0,'cubic coefficient');
    b=M(1:end-1)/2;normal(b,M(1:end-1)~=0,'quadratic coefficient');
    weighted=(2*M(1:end-1)+M(2:end))/6;
    normal(weighted,weighted~=0,'weighted second derivative');
    correction=h.*weighted;normal(correction,weighted~=0,'linear correction');
    c=delta-correction;normal(c,c~=0,'linear coefficient');
    coefficients=[a,b,c,q(1:end-1)];
    normal(coefficients,coefficients~=0,'coefficient');
    % Check the original local coefficients, not only the linear equations.
    left=[q(1:end-1),c,2*b];
    right=[((a.*h+b).*h+c).*h+q(1:end-1),(3*a.*h+2*b).*h+c,6*a.*h+2*b];
    if any(~isfinite([left(:);right(:)])),unresolved('Endpoint arithmetic is nonfinite.');end
    if any(abs(right(:,1)-q(2:end))>1e-9)
        unresolved('Stored coefficients do not reliably interpolate samples.');
    end
    join_tol=[1e-9,1e-8,1e-7]+1e-10*max(abs([left;right]),[],1);
    if n>2&&any(any(abs(right(1:end-1,:)-left(2:end,:))>join_tol))
        unresolved('Stored coefficient joins are unreliable.');
    end
    if (~isempty(v0)&&abs(left(1,2)-double(v0))>1e-8)||...
            (~isempty(vf)&&abs(right(end,2)-double(vf))>1e-8)
        unresolved('Stored coefficients do not reliably satisfy endpoint velocity.');
    end
    pp=mkpp(t',coefficients);
end
function normal(value,nonzero,where)
    if any(~isfinite(value(:)))||any(nonzero(:)&(value(:)==0|abs(value(:))<realmin))
        unresolved(['Unrepresentable ',where,'.']);
    end
end
function unresolved(message)
    error('coneik_boundary_spline:NumericalUnresolved','%s',message);
end
