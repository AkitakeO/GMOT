function [Mhat,FNMhat,Weight]=PROJonBRSURF(Shxyz,Nhxyz,Bhxyz,Chxyz0,ChInform,MEG2HEAD)
% Project senosrs on the brains surface
% Discussed in the HEAD coordinate space
% Shxyz: 3 x 102 (Positin of the sensors)
% Nhxyz: 3 x 102 (Radial unit vector of the sensors)
% Bhxyz: 3 x n (Position of the brain mesh points)
% Chxyz0: 3 x 1 (Position of the center of the brain)
% ChInform: 9 x 306 (Information of the sensors)
% MEG2HEAD: 3 x 8 (Converting matrix from MEG to HEAD coortidate)
% Mhat: 4 x 102 (4: x +y + z + topography value)
% FNMhat: 4 x 918
% Weight: n x 918
Mhat=zeros(4,size(Shxyz,2)); % 4x102
for s=1:size(Shxyz,2)
    yy=[Bhxyz(1,:)-repmat(Shxyz(1,s),[1,size(Bhxyz,2)]);...
        Bhxyz(2,:)-repmat(Shxyz(2,s),[1,size(Bhxyz,2)]);...
        Bhxyz(3,:)-repmat(Shxyz(3,s),[1,size(Bhxyz,2)])];
    H=[Nhxyz(1,s);Nhxyz(2,s);Nhxyz(3,s)]; % 3x1
    ddhat=(H'*H)\H'*yy;
    D2=sum((yy-H*ddhat).^2,1);
    yr=[Chxyz0(1)-Shxyz(1,s); Chxyz0(2)-Shxyz(2,s); Chxyz0(3)-Shxyz(3,s)];
    drhat=(H'*H)\H'*yr;
    ModD2=D2.*(ddhat.^2<drhat^2);
    ModD2=ModD2+(ddhat.^2>=drhat^2).*max(D2); 
    [~,Ihatsor]=sort(ModD2);
    dhatsor=ddhat(Ihatsor(1:50));
    [dhat,I]=max(dhatsor);
    Ihat=Ihatsor(I);
    D2hat=ModD2(Ihat);
    if D2hat>1e-4
        if f.gmftoption(1)==1
            Shxyz(:,s)=-Nhxyz(:,s)*0.07+Shxyz(:,s);
        elseif f.gmftoption(2)==2
            Shxyz(3,s)=Shxyz(3,s)-10;
        else
            Shxyz(:,s)=Nhxyz(:,s)*dhat+Shxyz(:,s);          
        end
        Nhxyz(:,s)=(Shxyz(:,s)-Chxyz0(1:3))/norm(Shxyz(:,s)-Chxyz0(1:3));
        yy=[Bhxyz(1,:)-repmat(Shxyz(1,s),[1,size(Bhxyz,2)]);...
            Bhxyz(2,:)-repmat(Shxyz(2,s),[1,size(Bhxyz,2)]);...
            Bhxyz(3,:)-repmat(Shxyz(3,s),[1,size(Bhxyz,2)])];
        H=[Nhxyz(1,s);Nhxyz(2,s);Nhxyz(3,s)]; % 3x1
        ddhat=(H'*H)\H'*yy;
        D2=sum((yy-H*ddhat).^2,1);
        ModD2=D2.*(ddhat.^2<drhat^2);
        ModD2=ModD2+(ddhat.^2>=drhat^2).*max(D2);
        [~,Ihatsor]=sort(ModD2);
        dhatsor=ddhat(Ihatsor(1:15));
        [~,I]=max(dhatsor);
        Ihat=Ihatsor(I);
    end
    
    Mhat(1:3,s)=Shxyz(:,s)+Nhxyz(:,s)*ddhat(Ihat);
end

cent=ChInform(1:3,3:3:306); % 3x102. Positon of the sensors (MEG cooridate)
unix=ChInform(4:6,3:3:306); % 3x102. Unit vectors in the x direction (MEG coordinate)
uniy=ChInform(7:9,3:3:306); % 3x102. Unit vectors in the y direction (MEG coordinate)
FNpointsS=repelem(cent,1,9);
FX=repelem(unix,1,9);
modx=repmat(-0.0045: 0.0045: 0.0045,[3,3*102]);
FNX=FX.*modx;
FY=repelem(uniy,1,9);
mody=repmat(-0.0045: 0.0045: 0.0045,[3,1]);
mody=reshape(mody,[1,9]);
mody=repmat(mody,[3,102]);
FNY=FY.*mody;
FNpointsS=FNpointsS+FNX+FNY;
ConvSH=MEG2HEAD;
FNpointsH=ConvSH(:,1:3)*FNpointsS+ConvSH(:,4);
Shxyz=FNpointsH;
Ndxyz=ChInform(10:12,3:3:306);
Nhxyz=ConvSH(:,1:3)*Ndxyz;
Nhxyz=repelem(Nhxyz,1,9);
FNMhat=zeros(4,size(Shxyz,2));
for s=1:size(Shxyz,2)
    yy=[Bhxyz(1,:)-repmat(Shxyz(1,s),[1,size(Bhxyz,2)]);...
        Bhxyz(2,:)-repmat(Shxyz(2,s),[1,size(Bhxyz,2)]);...
        Bhxyz(3,:)-repmat(Shxyz(3,s),[1,size(Bhxyz,2)])]; 
    H=[Nhxyz(1,s);Nhxyz(2,s);Nhxyz(3,s)]; % 3x1
    ddhat=(H'*H)\H'*yy;
    D2=sum((yy-H*ddhat).^2,1);
    yr=[Chxyz0(1)-Shxyz(1,s); Chxyz0(2)-Shxyz(2,s); Chxyz0(3)-Shxyz(3,s)];
    drhat=(H'*H)\H'*yr;
    ModD2=D2.*(ddhat.^2<drhat^2);
    ModD2=ModD2+(ddhat.^2>=drhat^2).*max(D2);
    [~,Ihatsor]=sort(ModD2);
    dhatsor=ddhat(Ihatsor(1:15));
    [dhat,I]=max(dhatsor);
    Ihat=Ihatsor(I);
    D2hat=ModD2(Ihat);
    if D2hat>1e-4
        if f.gmftoption(1)==1
            Shxyz(:,s)=-Nhxyz(:,s)*0.07+Shxyz(:,s);
        elseif f.gmftoption(2)==2
            Shxyz(3,s)=Shxyz(3,s)-10;
        else
            Shxyz(:,s)=Nhxyz(:,s)*dhat+Shxyz(:,s);
        end
        Nhxyz(:,s)=(Shxyz(:,s)-Chxyz0(1:3))/norm(Shxyz(:,s)-Chxyz0(1:3));
        yy=[Bhxyz(1,:)-repmat(Shxyz(1,s),[1,size(Bhxyz,2)]);...
            Bhxyz(2,:)-repmat(Shxyz(2,s),[1,size(Bhxyz,2)]);...
            Bhxyz(3,:)-repmat(Shxyz(3,s),[1,size(Bhxyz,2)])]; 
        H=[Nhxyz(1,s);Nhxyz(2,s);Nhxyz(3,s)];
        ddhat=(H'*H)\H'*yy;
        D2=sum((yy-H*ddhat).^2,1);
        ModD2=D2.*(ddhat.^2<drhat^2);
        ModD2=ModD2+(ddhat.^2>=drhat^2).*max(D2); 
        [~,Ihatsor]=sort(ModD2);
        dhatsor=ddhat(Ihatsor(1:15));
        [~,I]=max(dhatsor);
        Ihat=Ihatsor(I);
    end   
    FNMhat(1:3,s)=Shxyz(:,s)+Nhxyz(:,s)*ddhat(Ihat);
end
Weight=ones(size(Bhxyz,2),size(R.FNMhat,2)); % 2496x49*102
for s=1:size(R.FNMhat,2) % 102
    dist2=sum((Bhxyz-repmat(R.FNMhat(1:3,s),[1,size(Bhxyz,2)])).^(f.gmftoption(4)*2)); % 1x2496
    Weight(:,s)=1./dist2';
end