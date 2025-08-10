function ret=binom(n,m)
    ret=repmat(zeros('like',n)+zeros('like',m),[length(n) length(m)]);
    for i1=1:length(n)
        n_=n(i1);
        for i2=1:length(m)
            m_=m(i2);
            if ~isa(n_,'sym')&&n_<m_
                ret_=0;
            else
                ret_=1+0*n_;
                for ii=0:m_-1
                    ret_=(n_-ii)/(m_-ii)*ret_;
                end
            end
            ret(i1,i2)=ret_;
        end
    end

end