function ret=binom(n,m)
    ret=(n-n)+(m-m);
    n=repmat(n,size(ret)./size(n));
    m=repmat(m,size(ret)./size(m));
    ret=arrayfun(@binom_,n,m);
    function ret_=binom_(n_,m_)
        if ~isa(n_,'sym')&&n_<m_
            ret_=0;
        else
            ret_=1+0*n_;
            for ii=0:m_-1
                ret_=(n_-ii)/(m_-ii)*ret_;
            end
        end
    end

end