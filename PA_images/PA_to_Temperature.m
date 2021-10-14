% Conversion to temperature again:
T_med = zeros(3,col,t_size);
T_med(3,:,1) = 37;
for j=1:t_size
    T_med(1:2,:,j) = PA_0(1:2,:);
end
% For to complete T_m matrix in time: T_m = T_m0*PA/PA_0 con T_0=37C
for k=1:col
    T_med(3,k,2:t_size)=37*PA_t(3,k,2:t_size)/PA_0(3,k);
end

% Conversion to temperature again with movement:
T_medm = zeros(3,a*b,t_size);
T_medm(3,:,1) = 37;
for j=1:t_size
    T_medm(1:2,:,j) = PA_window(1:2,:);
end
% For to complete T_m matrix in time: T_m = T_m0*PA/PA_0 con T_0=37C
for k=1:a*b
    T_medm(3,k,2:t_size)=37*PA_tm(3,k,2:t_size)/PA_window(3,k);
end

%% Signal Processing:
% Filtering data:
T_filter = zeros(3,a*b,t_size);
T_filter(1:2,:,:) = T_medm(1:2,:,:);
T_filter(3,:,1) = 37;
PA_filter = zeros(3,a*b,t_size);
PA_filter(1:2,:,:) = T_medm(1:2,:,:);
%[coefb,coefa] = butter(6,1/(4/2));
filt = tf(1, [2/4 1]); %Filter (first order TF with a settling time around 2.5 sec. This means that fast changes won't affect the system
%for example, noise; but movement changes will.
filt = ss(filt);
for k=1:a*b
    %Tinfo = T_medm(3,k,1:t_size);
    %Tinfo = reshape(Tinfo,1,t_size);
    PAinfo = PA_tm(3,k,1:t_size);
    PAinfo = reshape(PAinfo,1,t_size);
    %T_filter(3,k,1:t_size) = filtfilt(coefb,coefa,Tinfo);
    %T_filter(3,k,1:t_size) = lsim(filt,Tinfo,tlist,37);
    PA_filter(3,k,1:t_size) = lsim(filt,PAinfo,tlist,PA_tm(3,k,1));
end
%Calculo de T_filter (esto es util solo para ploteo porque todo el
%algoritmo procesa PA y va calculando T con eso)
for k=1:a*b
    T_filter(3,k,2:t_size)=37*PA_filter(3,k,2:t_size)/PA_window(3,k);
end

[~,nodo]= min((PA_window(1,:) - 0.30).^2 + (PA_window(2,:) - 0.15).^2);
%cat(2, error_TD, 100-TD(end));
Tdata = zeros(3,a*b,t_size+2);
Tdata(1:2,:,:) = T_medm(1:2,:,:);
Tdata(3,:,1:2) = 37;
T_delay = 2; %Tdata is going to have a delay of 2n (0.5 s).
%Tdata = reshape(T_filter(3,nodo,1:4),1,4);
%avg = mov_ave(Tdata,0,1,2);
condnewgr = false; %Condition that states if the current state doesn't follow the average of the last data.
for tm=3:4
    if(tm==3)
        Tdata(3,:,tm) = T_filter(3,:,tm-2);
        avg = mov_ave(Tdata(3,:,:),a*b,0,2,tm);
        Tdata(3,:,tm) = avg; 
    elseif(tm==4)
        %Define the Reference node:
        [~,ref_node]=max(T_filter(3,:,tm-2));
        T_cur_ref = reshape(T_filter(3,ref_node,tm-2),1,1);
        ref_avg = avg(1,ref_node);
        T_cur_avg = (ref_avg+T_cur_ref)/2;
        if(T_cur_avg<=1+ref_avg)&&(T_cur_avg>=ref_avg-1)%Follows the expected results
            Tdata(3,:,tm) = T_filter(3,:,tm-2);
            avg = mov_ave(Tdata(3,:,:),a*b,0,2,tm);
            Tdata(3,:,tm) = avg; 
        else
            %Review if it is noise or not
            T_next_avg = mov_ave_uni(T_filter(3,ref_node,tm-1:tm),0,1,2);
            if(T_next_avg<=1.5+ref_avg)&&(T_next_avg>=ref_avg-1.5)%It is noise
                avg_next = mov_ave(T_filter(3,:,:),a*b,0,1,tm);
                Tdata(3,:,tm) = (T_next_avg+avg)/2;
                avg = mov_ave(Tdata(3,:,:),a*b,0,2,tm);
                Tdata(3,:,tm) = avg;
            else %It's not noise
                if(abs(T_next_avg-T_cur_ref)<abs(T_cur_avg-T_cur_ref))%Current T has more relationship with the next T values than with the past values
                    %Sigue buscar el nuevo pixel correspondiente
                    near_val = zero(2,8);
                    cur_x = reshape(T_filter(1,ref_node,tm-3),1,1);
                    cur_y = reshape(T_filter(2,ref_node,tm-3),1,1);
                    yval = cur_y-0.5;
                    contador = 1;
                    while yval <= cur_y+0.5
                        xval = cur_x-0.5;
                        while xval <= cur_x+0.5
                            if(yval~=cur_y)||(xval~=cur_x)
                                %Search for surrounding nodes and their
                                %average PA signal in the next two 
                                [~,nearNode]= min((PA_window(1,:) - xval).^2 + (PA_window(2,:) - yval).^2);
                                PA_ave_node = reshape(PA_filter(3,nearNode,tm-3),1,1); %mov_ave_uni(PA_filter(3,nearNode,tm-4:tm-3),0,1,2); %especificar caso de tm=5
                                near_val(:,contador) = [nearNode;PA_ave_node];
                            end
                            xval = xval+0.5;
                            contador = contador + 1;
                        end
                    end
                    
                else %T belongs to previous group
                    T_cur_ref = interp1([0.25 0.5],Tdata(tm-2:tm-1),[0.75],'linear','extrap');
                    Tdata = cat(2, Tdata(tm-3:tm-1),T_cur_ref);
                    avg = mov_ave(Tdata,0,3,tm);
                end
            end
        end
        
    end
    
end


% Conversion to temperature again with filter:
[X,Y]=meshgrid(0:pixel_size:1,-1:pixel_size:1);
[sizey,sizex]=size(X);
T_presmooth = zeros(sizey,sizex,t_size);
for i=1:sizex
    actual_x = X(1,i);
    for j=1:sizey
        actual_y = Y(j,1);
        [~,nodo]= min((T_medm(1,:,1) - actual_x).^2 + (T_medm(2,:,1) - actual_y).^2);
        T_presmooth(j,i,:)=T_medm(3,nodo,:);
    end
end

%Filter:
h = fspecial('average'); %3x3 filter
%PA signal smooth:
T_smooth = zeros(sizey,sizex,t_size);
for i=1:t_size
   T_smooth(:,:,i) = imfilter(T_presmooth(:,:,i),h,'replicate'); 
end
%Reshaping:
T_smooth_reshaped = zeros(3,a*b,t_size);
T_smooth_reshaped(1:2,:,:) = PA_tm(1:2,:,:);
for c=1:a*b
    Y_t = transpose(Y);
    [~,row]=min((Y_t(1,:)-T_smooth_reshaped(2,c,1)).^2);
    [~,colm]=min((X(1,:)-T_smooth_reshaped(1,c,1)).^2);
    T_smooth_reshaped(3,c,:)=T_smooth(row,colm,:);
end



%Plot of same node for original and reconstructed Temperature.
% [~,nodo]= min((PA_0(1,:) - 0.30).^2 + (PA_0(2,:) - 0.15).^2);
% [~,nid] = getClosestNode(mesh.Nodes, 0.30,0.15);
% figure;
% plot(tlist, T(nid,:),'.');
% hold on;
% plotT =reshape(T_med(3,nodo,:),1,t_size);
% plot(tlist,plotT,'--');
% [~,nodo]= min((PA_window(1,:) - 0.30).^2 + (PA_window(2,:) - 0.15).^2);
% plotTm =reshape(T_medm(3,nodo,:),1,t_size);
% plot(tlist,plotTm);
% plotTF =reshape(T_filter(3,nodo,:),1,t_size);
% plot(tlist,plotTF);
% legend({'Original','Reconstructed','Movement','Filtered'});
% title 'Temperature reconstructed and original at (0.30,0.15)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
% hold off;
% 
% %% Plot PA curves
% [~,nodo]= min((PA_0(1,:) - 0.30).^2 + (PA_0(2,:) - 0.15).^2);
% figure;
% plotPA_t =reshape(PA_t(3,nodo,:),1,t_size);
% plot(tlist, plotPA_t,'.');
% hold on;
% plotPA_tn =reshape(PA_tn(3,nodo,:),1,t_size);
% plot(tlist,plotPA_tn,'--');
% [~,nodo]= min((PA_window(1,:) - 0.30).^2 + (PA_window(2,:) - 0.15).^2);
% plotPA_tm =reshape(PA_tm(3,nodo,:),1,t_size);
% plot(tlist,plotPA_tm);
% plotPAF =reshape(PA_filter(3,nodo,:),1,t_size);
% plot(tlist,plotPAF);
% legend({'Original','Reconstructed','Movement','Filtered'});
% title 'PA at (0.30,0.15)';
% xlabel 'Time, seconds'
% ylabel 'PA'
% hold off;


[a,b]=size(X);
T_XY = zeros(a,b,t_size);
%Z = 0;
for i=1:b
    actual_x = X(1,i);
    for j=1:a
        actual_y = Y(j,1);
        %[~,nodo]= min((T_med(1,:,800) - actual_x).^2 + (T_med(2,:,800) - actual_y).^2);
        [~,nodo]=getClosestNode(mesh.Nodes,actual_x,actual_y);
        %T_XY(j,i)= T_med(3,nodo,800); 
        T_XY(j,i,:)= T(nodo,:);
    end
end
figure;
surf(X,Y,T_XY);
colorbar;

%% Moving average function:
function av = mov_ave(arr,length_xy,M1,M2,n)
    av = zeros(1,length_xy);
    for l=1:length_xy
        av_element = 0;
        for k=-M1:M2
            av_element = av_element+arr(3,l,n-k);
        end
        av(1,l) = av_element/(M1+M2+1);
    end
end

function av_uni = mov_ave_uni(arr,M1,M2,n)
    av_uni = 0;
    for k=-M1:M2
        av_uni = av_element+arr(1,1,n-k);
    end
    av_uni = av_uni/(M1+M2+1);
end