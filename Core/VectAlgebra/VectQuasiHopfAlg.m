classdef VectQuasiHopfAlg<VectAlg
    % VectQuasiHopfAlg : Quasi-Hopf 代数を表すクラス
    %  VectAlg を継承し、Quasi-Hopf 代数に特有の構造を追加
    % quasi-coassociativity, antipode, associator, twister など
    % issue: left/right unitorは自明なものとして省略する
    properties
        bracketing (1,:)char
        dyck (:,1) logical
    end
    methods
        function setSCQuasi(obj,v)
            arguments
                obj
                v.associator SparseEx % associator Phi 
                v.associator_inv SparseEx % inverse of associator Phi^{-1}
                v.twister_F SparseEx % twister F-tensor
                v.alpha_S SparseEx % element α for antipode
                v.beta_S SparseEx % element β for antipode
                v.unitor_l SparseEx % left unit
                v.unitor_r SparseEx % right unit    
            end
            s=obj.spec;
            for fld=string(fieldnames(v))'
                s.SC{fld} = v.(fld);
            end
        end
        function verifyHopf(obj,isequal_,type)
            arguments
                obj
                isequal_ (1,1)
                type (1,:) string {mustBeMember(type,["bialg","hopf","twister","all"])}="all"
            end
            if type=="all"
                type=["bialg","hopf","twister"];
            end
            D=obj.dim;
            I=SparseEx(eye(D));
            mu=obj.getSC('prod');
            eta=obj.getSC('unit');
            Delta=obj.getSC('coprod');
            ep=obj.getSC('counit');
            S=obj.getSC('antipode');
            Si=obj.getSC('antipode_inv');
            Phi=obj.getSC('associator') ;
            Phi_inv=obj.getSC('associator_inv');
            if ismember("bialg",type)
                % size
                assert(isequal(mu.size,[D D D]),"invalid size:μ")
                assert(isequal(eta.size,[D]),"invalid size:η")
                assert(isequal(Delta.size,[D D D]),"invalid size:Δ")
                assert(isequal(ep.size,[D]),"invalid size:ε")
                
                % unitality
                % tmp=tensorprod(eta,mu,1,1,Num=1);
                tmp=calcTensorExpression('eta{1}mu{1,2,3}',[2,3]);
                assertT(tmp,I,"left unitality error")
                % tmp=tensorprod(eta,mu,1,2,Num=1);
                tmp=calcTensorExpression('eta{2}mu{1,2,3}',[1,3]);
                assertT(tmp,I,"right unitality error")
                % counitality
                % tmp=tensorprod(ep,Delta,1,2,Num=1);
                tmp=calcTensorExpression('ep{2}Delta{1,2,3}',[1,3]);
                assertT(tmp,I,"left counitality error")
                % tmp=tensorprod(ep,Delta,1,3,Num=1);
                tmp=calcTensorExpression('ep{3}Delta{1,2,3}',[1,2]);
                assertT(tmp,I,"right counitality error")
                
                % associativity
                % tmp=tensorprod(mu,mu,3,1);
                tmp=calcTensorExpression('mu{1,2,3}mu{3,4,5}',[1,2,4,5]);
                % tmp2=permute(tensorprod(mu,mu,2,3),[1 3 4 2]);
                tmp2=calcTensorExpression('mu{1,2,3}mu{4,5,2}',[1,4,5,3]);
                assertT(tmp,tmp2,"associativity error")
                
            end

            
            % quasi-coassociativity
            % tmp=permute(tensorprod(Delta,Delta,2,1),[1 3 4 2]);
            tmp=calcTensorExpression('Delta{1,2,3}Delta{2,4,5}',[1,4,5,3]);
            % tmp2=tensorprod(Delta,Delta,3,1);
            tmp2=calcTensorExpression('Delta{1,2,3}Delta{3,4,5}',[1,2,4,5]);
            assertT(tmp,tmp2,"coassociativity error")

            if ismember("twister",type)
                F=obj.getSC('twister_F');
                assert(isequal(F.size,[D D]),"invalid size:F")

                % twist property
                % Δ_F = F Δ F^{-1}
                % tmp=tensorprod(F,Delta,2,1);
                tmp=calcTensorExpression('F{1,2}Delta{2,3,4}',[1,3,4]);
                % tmp=tensorprod(tmp,Fi,2,2);
                tmp=calcTensorExpression('F{1,2}Delta{2,3,4}F{4,5}',[1,3,5]);
                % tmp2=Delta_F;
                tmp2=obj.getSC('coprod_twisted');
                assertT(tmp,tmp2,"twist coproduct error")
                % μ_F = F^{-1} μ F
                % tmp=tensorprod(mu,F,3,1);
                tmp=calcTensorExpression('mu{1,2,3}F{3,4}',[1,2,4]);
                % tmp=tensorprod(Fi,tmp,2,2);
                tmp=calcTensorExpression('F{1,2}mu{2,3,4}F{4,5}',[1,3,5]);
                % tmp2=mu_F;
                tmp2=obj.getSC('prod_twisted');
                assertT(tmp,tmp2,"twist product error")
            end
            if ismember("hopf",type)

            
                alpha=obj.getSC('alpha_S');
                beta=obj.getSC('beta_S');
                
                assert(isequal(S.size,[D D]),"invalid size:S")
                assert(isequal(Si.size,[D D]),"invalid size:Si")
                assert(isequal(alpha.size,[D]),"invalid size:α")
                assert(isequal(beta.size,[D]),"invalid size:β")
            end


            
            
            
            % bialgebra property
            tmp=calcTensorExpression( ...
                'Delta{1,3,4}Delta{2,5,6}mu{3,5,7}mu{4,6,8}',[1,2,7,8]);
            tmp2=calcTensorExpression( ...
                'mu{1,2,3}Delta{3,7,8}',[1,2,7,8]);
            assertT(tmp,tmp2,"bialgebra property error")
            
            
            % antipode property
            % μ ∘ (S ⊗ id) ∘ Δ = η ∘ ε  and  μ ∘ (id ⊗ S) ∘ Δ = η ∘ ε
            % tmp2=tensorprod(ep,eta,2,2,Num=2);
            tmp2=calcTensorExpression('ep{1}eta{2}',[1,2]);
            
            % tmp=tensorprod(S,Delta,2,2);
            % tmp=tensorprod(tmp,mu,[1,3],[1,2]);
            tmp=calcTensorExpression('S{2,3}Delta{1,3,4}mu{2,4,5}',[1,5]);
            assertT(tmp,tmp2,"antipode left inverse error")
            
            % tmp=tensorprod(Delta,S,3,2);
            % tmp=tensorprod(tmp,mu,[2 3],[1 2]);
            tmp=calcTensorExpression('Delta{1,2,3}S{4,3}mu{2,4,5}',[1,5]);
            assertT(tmp,tmp2,"antipode right inverse error")

            disp("Confirmed to be a Hopf algebra")
            function assertT(val1,val2,msg)
                diff=C(val1-val2);
                cond=isequal_(diff.val,0);
                if ~cond
                    error("VectAlg:verifyHopf",msg)
                end
            end
            
            % === quasi-Hopf 代数の標準公理（参考用コメント）=================
            % ここで扱っている構造:
            %  (H, μ, η, Δ, ε, Φ, S, α, β)
            %
            % 1) 弱コ多重性 (quasi-coassociativity):
            %    (Δ ⊗ id)Δ(h) · Φ  =  Φ · (id ⊗ Δ)Δ(h)
            %    （コード内では Φ による歪みを直接検証していない箇所あり）
            %
            % 2) 3-コサイクル条件 (pentagon):
            %    (1⊗Φ)(id⊗Δ⊗id)(Φ)(Φ⊗1) = (id⊗id⊗Δ)(Φ)(Δ⊗id⊗id)(Φ)
            %
            % 3) counit:
            %    (ε ⊗ id)Δ = id = (id ⊗ ε)Δ
            %
            % 4) S, α, β の公理（S は必ずしも反代数射でなくても良い）:
            %    S(h_1) α h_2 = ε(h) α
            %    h_1 β S(h_2) = ε(h) β
            %
            % 5) Φ = Σ X^1 ⊗ X^2 ⊗ X^3,  Φ^{-1} = Σ x^1 ⊗ x^2 ⊗ x^3 を用いて:
            %    Σ X^1 β S(X^2) α X^3 = 1
            %    Σ S(x^1) α x^2 β S(x^3) = 1
            %
            % 6) (H-Module 圏の単位制約を明示的に線形写像として格納したい場合):
            %    unit_l, unit_r を
            %       l_V : k ⊗ V → V,   r_V : V ⊗ k → V
            %    の行列表現（基底固定後）に対応させることができる。
            %    そのとき三角恒等式は associator a_{U,V,W}(=Φ 作用) を使って:
            %       (id_V ⊗ l_W) ∘ a_{V,k,W} = r_V ⊗ id_W
            %       (l_V ⊗ id_W) = (id_V ⊗ r_W) ∘ a_{V,W,k}
            %    etc.（厳密形は圏の取り扱い方針に依存）
            %
            % 7) 本クラスの unit_l / unit_r を「別個の単位射データ」として保持する場合、
            %    通常 Hopf では冗長だが、quasi-Monoidal な補正や基底変換を明示する目的で使う設計もあり得る。
            %
            % 8) もし unit_l / unit_r が不要（η だけで十分）なら SC セット時に省略し、
            %    依存する検証ロジックを追加しない方が一貫する。
            %
            % ※必要なら今後:
            %    - Φ pentagon の数値検証
            %    - α, β 公理の数値検証
            %    - unit_l / unit_r を用いた「モジュール圏的三角恒等式」の実装
            %   を追加可能。
            % =============================================================
        end
    end
end