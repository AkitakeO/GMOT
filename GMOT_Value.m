function Values=GMOT_Value(GTdata,GLdata,tstart,sampfreq)
% GTdata n x 102 (Data of transverse planar gradiometers)
% GLdata n x 102 (Data of longitudinal planar gradiometers)
% tstart: Time to start
% sampfreq: Sampling frequency = 1,000 Hz
p=8; % Getting 8 values
DeSMagni=zeros(p,size(GTdata,2)); % p*102
ThSMagni=zeros(p,size(GTdata,2)); % p*102
AlSMagni=zeros(p,size(GTdata,2)); % p*102
BeSMagni=zeros(p,size(GTdata,2)); % p*102
GaSMagni=zeros(p,size(GTdata,2)); % p*102
GHSMagni=zeros(p,size(GTdata,2)); % p*102
LHSMagni=zeros(p,size(GTdata,2)); % p*102
MHSMagni=zeros(p,size(GTdata,2)); % p*102
HHSMagni=zeros(p,size(GTdata,2)); % p*102
for m=1:p 
    tsamp=(round(((m-1)*0.5+tstart)*sampfreq)+1:round(((m-1)*0.5+tstart)*sampfreq)+2^10)'; % 1024x1
    SampGTdata=GTdata(tsamp,:);
    SampGLdata=GLdata(tsamp,:);
        SampGTdata=SVD_filter(SampGTdata);
        SampGLdata=SVD_filter(SampGLdata);    
    SMagni=zeros(2^9+1,size(R.GTdata,2)); 
    for n=1:102 
        magniT=spectrogram(SampGTdata(:,n),2^10,2^9); 
        magniT=abs(magniT).^2/2^10;
        SMagniT=mean(magniT,2);
        magniL=spectrogram(SampGLdata(:,n),2^10,2^9);
        magniL=abs(magniL).^2/2^10;
        SMagniL=mean(magniL,2);
        SMagni(:,n)=sqrt(SMagniT.^2+SMagniL.^2);
    end
    ghnum=[46:57 63:70];
    lhnum=71:100;
    mhnum=[101:117 123:177 183:200];
    hhnum=[201:217 223:237 243:297 303:330];
    
    DeSMagni(m,:)=mean(SMagni(round((0.5:3)/(sampfreq*0.5)*512)+1,:)); % 1x102<-3x102 (2-4th sample)
    ThSMagni(m,:)=mean(SMagni(round((4:7)/(sampfreq*0.5)*512)+1,:)); % 1x102<-4x102 (5-8th sample)
    AlSMagni(m,:)=mean(SMagni(round((8:13)/(sampfreq*0.5)*512)+1,:)); % 1x102<-6x102 (9-14th sample)
    BeSMagni(m,:)=mean(SMagni(round((14:30)/(sampfreq*0.5)*512)+1,:)); % 1x102<-18x102 (15-32nd sample)
    GaSMagni(m,:)=mean(SMagni(round((26:45)/(sampfreq*0.5)*512)+1,:)); % 1x102<-20x102 (28-47th sample)
    GHSMagni(m,:)=mean(SMagni(round((ghnum)/(sampfreq*0.5)*512)+1,:)); % 1x102<-20x102 (28-47th sample
    LHSMagni(m,:)=mean(SMagni(round((lhnum)/(sampfreq*0.5)*512)+1,:)); % 1x102<-20x102 (28-47th sample
    MHSMagni(m,:)=mean(SMagni(round((mhnum)/(sampfreq*0.5)*512)+1,:)); % 1x102<-20x102 (28-47th sample
    HHSMagni(m,:)=mean(SMagni(round((hhnum)/(sampfreq*0.5)*512)+1,:)); % 1x102<-20x102 (28-47th sample
end
Values=zeros([size(DeSMagni),9]);
Values(:,:,1)=DeSMagni;
Values(:,:,2)=ThSMagni;
Values(:,:,3)=AlSMagni;
Values(:,:,4)=BeSMagni;
Values(:,:,5)=GaSMagni;
Values(:,:,6)=GHSMagni;
Values(:,:,7)=LHSMagni;
Values(:,:,8)=MhSMagni;
Values(:,:,9)=HHSMagni;
end