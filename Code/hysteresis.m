function output_hyst=hysteresis(in_files, nor, diagram, file_el)
% 
% Determines hysteresis function from observation series at different
% stations from several days. Observation series should last for at least
% 30 min (cca).
%
% output_hyst=hysteresis(in_files, output, nor, diagram, file_el)
%
% in_files...  cell array with names and paths to input files with gravity
%              readings. One file = one measuring day.
% nor... numer of readings for each occupation to include in calculations.
%        If set to 0, the algorithm calculates with real number of readings
%        for each occupation (all readings).
% diagram = [0 1] -> Not plots or plots, respectively, the diagram with
%                    observation series for all occupations with adjusted
%                    hysteresis function.
% file_el = [0 1]  0 -> Doesn't creates data files with eliminated
%                       hysteresis effect.
%                  1 -> Creates data files with eliminated hysteresis.
%                       Original Scintrex format.
%                  2 -> Creates data files with eliminated hysteresis.
%                       Extended format with more decimal places.

%--------------------------------------------------------------------------
nod=size(in_files,1); % Number of days

% Initialisation of variables
gm_const_C=cell(nod,1);
data=cell(nod,1);
index_t=cell(nod,1);
z_sd_t_C=cell(nod,1);
o_n=zeros(nod, 1);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Data import
% -------------------------------------------------------------------------

for day=1:nod
    % Import of gravimeter's constants from file heading. Stored in
    % variable:
    % gm_const [DATE; GREF; GCAL1; GCAL2; TEMPCO; DRIFT_value; DRIFT_start_time;
    %           DRIFT_start_date; X SENS; Y SENS; GMT_diff; CYC_time; SER.NO;
    %           OPERATOR; FI; LA]
    % Conversion to micro ms^-2 is applied!
    gm_const=zeros(14,1);
    aa=fopen(in_files{day},'r');
    POSTAVKE = textscan(aa, 'Cycle Time: %f Ser No: %f\n', 1, 'headerlines', 3);
    gm_const(12)=POSTAVKE{1,1};
    gm_const(13)=POSTAVKE{1,2}; % multifunkcijska pomocna varijabla
    POSTAVKE = textscan(aa, 'Line: %*d Grid: %*d Job: %*d Date: %s Operator: %d\n\n', 1);
    p_str=cell2mat(POSTAVKE{1,1}); % multifunkcijska pomocna varijabla
    gm_const(1)=datenum([p_str(4:8) '/' p_str(1:2)]);
    gm_const(14)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'GREF.: %f mGals Tilt x sensit.: %f\n', 1);
    gm_const(2)=POSTAVKE{1,1}*10;
    gm_const(9)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'GCAL.1: %f Tilt y sensit.: %f\n', 1);
    gm_const(3)=POSTAVKE{1,1}*10;
    gm_const(10)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'GCAL.2: %f Deg.Latitude: %f\n', 1);
    gm_const(4)=POSTAVKE{1,1}*100;
    gm_const(15)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'TEMPCO.: %f mGal/mK Deg.Longitude: %f\n', 1);
    gm_const(5)=POSTAVKE{1,1}*10;
    gm_const(16)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'Drift const.: %f GMT Difference: %fhr\n', 1);
    gm_const(6)=POSTAVKE{1,1}*10;
    gm_const(11)=POSTAVKE{1,2};
    POSTAVKE = textscan(aa, 'Drift Correction Start  Time: %s Cal.after x samples: %d\n', 1);
    p_str=cell2mat(POSTAVKE{1,1});
    gm_const(7)=datenum(['0-jan-0000,' p_str]);
    POSTAVKE = textscan(aa, 'Date: %s', 1);
    p_str=cell2mat(POSTAVKE{1,1});
    gm_const(8)=datenum([p_str(4:8) '/' p_str(1:2)]);
    fclose(aa);
    gm_const_C{day}=gm_const;
    
    % Import of gravity readings
    % in_files{day} % Useful to find an error in data files
    data_cell=unos_mix(in_files{day},10,15,[2,10]);
        
    % Conversion from textual to numerical variables.
    g=strrep(data_cell{2},'*','');
    g=str2num(cell2mat(g));
    time=cell2mat(data_cell{10});
    time_s=(str2num(time(:,1:2))*60+str2num(time(:,4:5)))*60+str2num(time(:,7:8));  % vrijeme u sekundama
    data{day}=[data_cell{1} g data_cell{3:9} time_s/(24*60*60)]; % vrijeme je u danima!!
        
    % Sorting of readings according to time.
    [time_s,index_t{day}]=sort(time_s);
    data{day}=data{day}(index_t{day},:);
    
    % Determination of number of occupations.
    d_time=[time_s(2:size(time_s,1))-time_s(1:(size(time_s,1)-1)); time_s(size(time_s,1))-time_s(1)];
    %i=find(d_time>gm_const(12)*3); % Option
    i=find(d_time>600); % 10 min - Minimal time difference between occupations
    o_n(day)=size(i,1);
    
    % Initialisation
    z_sd_t_DAY=cell(o_n(day),1);
    
    % Collecting data from all days into a cell array
    i=[0;i];
    for j=1:o_n(day)
        
        %----------------------------------------------------------------------
        % Conversion is applied form mGal to micro ms-2 for gravity!
        % Conversion is applied form day to min for time!
        %----------------------------------------------------------------------

        z_sd_t_DAY{j}=[data{day}((i(j)+1):i(j+1),2)*10 data{day}((i(j)+1):i(j+1),3)*10 data{day}((i(j)+1):i(j+1),10)*24*60]; % time in min!
        z_sd_t_DAY{j}(:,3)= z_sd_t_DAY{j}(:,3)- z_sd_t_DAY{j}(1,3); % time is conting from the first reading

        if nor && nor<size(z_sd_t_DAY{j},1)
            z_sd_t_DAY{j}=z_sd_t_DAY{j}(1:nor,:);
        end

    end
    
    z_sd_t_C{day}=z_sd_t_DAY;
    
end

%--------------------------------------------------------------------------
% Preparation for adjustment.
%--------------------------------------------------------------------------

% y =a(k) + b * exp(c*(t-d(k)))

% a(k) is not determined in adjustment. The minimal value of z is taken and
% in each iteration the value is corrected until the trend in residuals can
% be detected.

% Conversion from cell array to array; forming the reference array:
% ref[day occupation]
z_sd_t_C1=cell(size(z_sd_t_C));
ref=[];

u=2;
for day=1:nod
    for st=1:o_n(day)
        br_z=size(z_sd_t_C{day}{st},1);
        u=u+1;
        ref=[ref; repmat([day st], br_z, 1)];
    end
    z_sd_t_C1{day}=cell2mat(z_sd_t_C{day});
end
z_sd_t=cell2mat(z_sd_t_C1);

n=size(z_sd_t,1);

index=find(ref(2:end,2)-ref(1:(end-1),2));
index=[0; index; size(ref,1)]; % indexes of end of data for specific occupation

%--------------------------------------------------------------------------
% Initialisation of variables and values for iteration 1.

% A, x0, ln (v=Ax-l, vtPv->min)
A=zeros(n,u);
A(:,1)=ones(n,1);
L=zeros(n,1);
ln=zeros(n,1);
p=zeros(n,1);

a0=zeros(size(index,1)-1,1);
Da0=zeros(size(index,1)-1,4);
x0=zeros(u,1);
x_def=zeros(u,1);

% -------------------------------------------------------------------------
% Approximate values of coefficients of exponential function. For
% increasing hysteresis or other types of gravimeters it could be necessary
% to adjust it.
b0=0.50;
c0=-0.05;

x_def(1:2)=[log(b0); c0];

% OFFSET for first iteration
for i=1:(size(index,1)-1)
    Da0(i,1)=min(z_sd_t((index(i)+1):index(i+1),1)); % minimal reading in observation series
end

%--------------------------------------------------------------------------
% Iterative adjustment - determination of final value of ak
%--------------------------------------------------------------------------


max_ia=[200 400 400];
step=[0.01 0.001 0.0001];

for k=1:3
    iteracija_a=0;
    while sum(round(abs(Da0(:,1)+Da0(:,2)+Da0(:,3)+Da0(:,4))*10^4)) && iteracija_a < max_ia(k)
        iteracija_a= iteracija_a+1;

        % -----------------------------------------------------------------
        % Initial adjustment x0_1

        index_poz=zeros(n,1);

        a0=a0+Da0(:,1);
        Da0(:,2:4)=Da0(:,1:3);
        Da0(:,1)=zeros(size(index,1)-1,1);

        x0=x_def;

        for i=1:(size(index,1)-1)
            index_poz((index(i)+1):index(i+1))=(z_sd_t((index(i)+1):index(i+1),1)-a0(i))>0;
            br_z=index(i+1)-(index(i)+1)+1;
            L((index(i)+1):index(i+1))=log(z_sd_t((index(i)+1):index(i+1),1)-a0(i));
            ln((index(i)+1):index(i+1))=x0(1)+x0(2)*(z_sd_t((index(i)+1):index(i+1),3)-x0(2+i))-L((index(i)+1):index(i+1));
            p((index(i)+1):index(i+1))=(exp(x0(1))*exp(x0(2)*( z_sd_t((index(i)+1):index(i+1),3) - x0(2+i)))).^2; % p = [b0 exp(c0 * (t - d0))]^2
            A((index(i)+1):index(i+1),2)=z_sd_t((index(i)+1):index(i+1),3)-x0(2+i);
            A((index(i)+1):index(i+1),2+i)=-1*ones(br_z,1)*x0(2);
        end

        % Elimination of non-positive values for L (ln(y) is not defined for
        % y<=0)
        index_poz=logical(index_poz);
        ref_r=ref(index_poz,:);
        index_r=find(ref_r(2:end,2)-ref_r(1:(end-1),2));
        index_r=[0; index_r; size(ref_r,1)]; % contains the ends of data for specific occupation

        L_r=L(index_poz);
        ln_r=ln(index_poz);
        p_r=p(index_poz);
        A_r=A(index_poz,:);

        %p_r=ones(size(p_r)); % For stochastic model P=I
        P=diag(p_r);

        x=-1*pinv(A_r'*P*A_r)*(A_r'*P*ln_r);
        v=A_r*x+ln_r;
        x_def=x0+x;

        % -----------------------------------------------------------------
        % Iterative adjustment - substitution with more accurate x0 because
        % of linearization of observation equations

        iteracijaX=0;
        while sum(abs(round(x*10^5)))
            iteracijaX=iteracijaX+1;
            x0=x_def;
            for i=1:(size(index,1)-1)
                br_z=index(i+1)-(index(i)+1)+1;
                ln((index(i)+1):index(i+1))=x0(1)+x0(2)*(z_sd_t((index(i)+1):index(i+1),3)-x0(2+i))-L((index(i)+1):index(i+1));
                p((index(i)+1):index(i+1))=(exp(x0(1))*exp(x0(2)*( z_sd_t((index(i)+1):index(i+1),3) - x0(2+i)))).^2; % p = [b0 exp(c0 * (t - d0))]^2
                A((index(i)+1):index(i+1),2)=z_sd_t((index(i)+1):index(i+1),3)-x0(2+i);
                A((index(i)+1):index(i+1),2+i)=-1*ones(br_z,1)*x0(2);
            end

            % Elimination of non-positive values for L (ln(y) is not
            % defined for y<=0)
            ln_r=ln(index_poz);
            A_r=A(index_poz,:);
            p_r=p(index_poz);
            %p_r=ones(size(p_r)); % For stochastic model P=I
            P=diag(p_r);
            x=-1*pinv(A_r'*P*A_r)*(A_r'*P*ln_r);
            v=A_r*x+ln_r;
            x_def=x0+x;
        end

        % Transformation to original function
        b=exp(x_def(1));
        c=x_def(2);

        % Determination of corrections for a(k)
        for i=1:(size(index,1)-1)
            % All measurements (negative and zero values also) are included
            % v = function - reading = b*exp((c*(t-d)) - z + a0
            v_org_st=b*exp(c*(z_sd_t((index(i)+1):index(i+1),3)-x_def(2+i)))-z_sd_t((index(i)+1):index(i+1),1)+a0(i);

            R = corrcoef(z_sd_t((index(i)+1):index(i+1),3),v_org_st);
            r=R(1,2);
            if abs(r)>0.01
                Da0(i,1)=-step(k)*r/abs(r);
            end

        end
    end
    
    Da0(:,2:4)=zeros(size(Da0(:,2:4)));
end

% Accuracy assessment

s0=sqrt((v'*P*v)/(size(A_r,1)-size(A_r,2)+1));
Qxx=pinv(A_r'*P*A_r);
sx=s0*sqrt(diag(Qxx));

% Transformation to original function
sb=b*sx(1);
sc=sx(2);

%--------------------------------------------------------------------------
% Diagram of adjusted hysteresis function with all observation series
%--------------------------------------------------------------------------
if diagram
    for i=1:(size(index,1)-1)
        plot((z_sd_t((index(i)+1):index(i+1),3)-x_def(2+i)),(z_sd_t((index(i)+1):index(i+1),1)-a0(i)))
        set(gca, 'NextPlot', 'Add');
    end
    t=-10:(1/60):70;
    f=b*exp(c*t);
    plot(t,f,'r');
    set(gca, 'NextPlot', 'Add');
    set(gca,'YLim',[-0.1 1.5])
    set(gca,'XLim',[-10 70])
    text(10,1.2,['y = (' num2str(b,'%5.3f') ' \pm ' num2str(sb,'%5.3f') ') * exp(' num2str(c,'%6.3f') ' \pm ' num2str(sc,'%5.3f') ' t)']);
    plot([-20,100],[0,0],'g');
    set(gcf, 'NextPlot', 'New');
end

%--------------------------------------------------------------------------
% Packing the data into an output variable.
%--------------------------------------------------------------------------
output_hyst=[x_def sx];
output_hyst(1:2, :)=[b sb; c sc];

%--------------------------------------------------------------------------
% Creation of data files with eliminated hysteresis effect.
%--------------------------------------------------------------------------
if file_el
    k=0;
    for day=1:nod
        h=0;
        for st=1:o_n(day)
            nor_st=size(z_sd_t_C{day}{st},1);
            data{day}(h+(1:nor_st),2)=(z_sd_t_C{day}{st}(:,1)-b*exp(c*(z_sd_t_C{day}{st}(:,3)-x_def(3+k))))/10;
            k=k+1;
            h=h+nor_st;
        end
        
        % sorting according to initial order
        data{day}(index_t{day},:)=data{day};

        %--------------------------------------------------------------------------
        % Creation and printout of new data file with eliminated hysteresis effect
        %--------------------------------------------------------------------------
        a=fopen([in_files{day}(1:end-4) '_ELH' in_files{day}(end-3:end)],'wt');

        % Printout of file header.

        fprintf(a,'\n-------------------------------------------------------------------------------\n');
        fprintf(a,'SCINTREX V5.2       AUTOGRAV / Cycling Mode          R5.31\n');
        fprintf(a,'Cycle Time:    %4u                                          Ser No:    %6u.\n',gm_const_C{day}(12),gm_const_C{day}(13));
        fprintf(a,'Line:   9999.  Grid:   9999.   Job:   9999.  Date: %s  Operator:      %2u.\n\n',datestr(gm_const_C{day}(1),25),gm_const_C{day}(14));
        fprintf(a,'GREF.:                 %8.3f mGals           Tilt x sensit.:           %5.1f\n',gm_const_C{day}(2)/10,gm_const_C{day}(9));
        fprintf(a,'GCAL.1:                %8.3f                 Tilt y sensit.:           %5.1f\n',gm_const_C{day}(3)/10,gm_const_C{day}(10));
        fprintf(a,'GCAL.2:                %8.3f                 Deg.Latitude:            %6.2f\n',gm_const_C{day}(4)/100,gm_const_C{day}(15));
        fprintf(a,'TEMPCO.:                %7.4f mGal/mK         Deg.Longitude:           %6.2f\n',gm_const_C{day}(5)/10,gm_const_C{day}(16));
        fprintf(a,'Drift const.:             %6.4f                GMT Difference:           %2d.hr\n',gm_const_C{day}(6)/10,gm_const_C{day}(11));
        fprintf(a,'Drift Correction Start  Time: %s          Cal.after x samples:         12\n',datestr(gm_const_C{day}(7),13));
        fprintf(a,'                        Date: %s          On-Line Tilt Corrected = "*"\n',datestr(gm_const_C{day}(8),25));
        fprintf(a,'-------------------------------------------------------------------------------\n');
        fprintf(a,'Station  Grav.     SD.   Tilt x  Tilt y    Temp.    E.T.C.  Dur  # Rej     Time\n');

        % Printout of readings 

        switch file_el
            case 1
                % Original format
                for i=1:size(data{day},1)
                    fprintf(a,'%6u. %8.3f* %5.3f   %4d.   %4d.    %4.2f   %6.3f   %4u  %4u  %s\n',data{day}(i,1:9),datestr(data{day}(i,10),13));
                end
            case 2
                for i=1:size(data{day},1)
                    % Extended format with extra decimal places
                    fprintf(a,'%5u. %9.4f* %5.3f   %4d.   %4d.   %5.2f  %7.4f   %4u  %4u  %s\n',data{day}(i,1:9),datestr(data{day}(i,10),13));
                end
        end
        fprintf(a,'\n');

        fclose(a);
    end
end
