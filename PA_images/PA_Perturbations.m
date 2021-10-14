t_move = [4];
condition = 0;
while(condition==0)
    t_new = randsample([4,5,6,7,8,9,10],1)+t_move(end);
    if(t_new>=300)
        condition = 1;
    else
        t_move = [t_move,t_new];
    end
end

for t=2:t_size
    % 1=up, 2=right, 3=down, 4=left
    %r = randsample([1 2 3 4],1);
    %First, we compare the corners of inner and outer squares:
    %Case 1: Top left corner. Node 41
    actualtime=tlist(t);
    [value,~] = min((t_move-actualtime).^2);
    if(value==0)
            paso = randsample([0.05 0.1 0.15 0.2],1);
            if (PA_tm(1,41,t-1)-paso<0) && (PA_tm(2,41,t-1)+paso>1.2)
                r = randsample([2 3],1);
                if r==2
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                end
            
            %Case 2: Buttom left corner. Node 1
            elseif (PA_tm(1,1,t-1)-paso<0) && (PA_tm(2,1,t-1)-paso<-1.2)
                r = randsample([1 2],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
                
            %Case 3: Top right corner. Node 861
            elseif (PA_tm(1,861,t-1)+paso>1.2) && (PA_tm(2,861,t-1)+paso>1.2)
                r = randsample([3 4],1);
                if r==3
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
            
            %Case 4: Buttom right corner. 821
            elseif (PA_tm(1,821,t-1)+paso>1.2) && (PA_tm(2,821,t-1)-paso<-1.2)
                r = randsample([1 4],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
                
            %We have to compare the 4 borders of the inner square vs the 4 borders
            %of the big square:
            %Case 5: moving the square more to the left would cause that it leaves
            %the tissue:
            elseif PA_tm(1,1,t-1)-paso<0
                %In this case, it could only be moved the other ways:
                r = randsample([1 2 3],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                elseif r==2
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                end
            
            %Case 6: moving the square more to the right would cause that it leaves
            %the tissue:
            elseif PA_tm(1,821,t-1)+paso>1.2
                r = randsample([1 3 4],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                elseif r==3
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
            
            %Case 7: moving the square upper would cause that it leaves
            %the tissue:
            elseif PA_tm(2,861,t-1)+paso>1.2
                r = randsample([2 3 4],1);
                if r==2
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                elseif r==3
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
            
            %Case 8: moving the square further down would cause that it leaves
            %the tissue:
            elseif PA_tm(2,821,t-1)-paso<-1.2
                r = randsample([1 2 4],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                elseif r==2
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
            
            else
                r = randsample([1 2 3 4],1);
                if r==1
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)+paso;
                elseif r==2
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)+paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                elseif r==3
                    PA_tm(1,:,t)=PA_tm(1,:,t-1);
                    PA_tm(2,:,t)=PA_tm(2,:,t-1)-paso;
                else
                    PA_tm(1,:,t)=PA_tm(1,:,t-1)-paso;
                    PA_tm(2,:,t)=PA_tm(2,:,t-1);
                end
            end 
    else
        PA_tm(1,:,t)=PA_tm(1,:,t-1);
        PA_tm(2,:,t)=PA_tm(2,:,t-1);
    end

    
    for k=1:a*b
        x1 = PA_tm(1,k,t);
        y1 = PA_tm(2,k,t);
        [~,node]= getClosestNode(PA_tn(:,:,1),x1,y1);
        PA_tm(:,k,t)=PA_tn(:,node,t);
    end
end