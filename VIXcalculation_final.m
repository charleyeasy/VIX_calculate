%% prepare data
% option_spx_hw2.csv : option prices for S&P 500 index
% rf_term_hw2.csv : risk-free rate for each term
% spx_price_hw2.csv : underlying asset (S&P500 index) price
%% S&P500 index options prices
% pre-import data (datastore)
dv = datastore('option_spx_final.csv');
% preview our data format
% preview(dv)
% select variables into our data
%%
dv.SelectedVariableNames = {'secid','date','exdate','cp_flag',...
  'strike_price','best_bid','best_offer','volume','ticker'};
% import data 
IVSPX = readall(dv);
% SPX data (SECID = 108105)
% strike_price = strike prcie/1000
IVSPX.strike_price = IVSPX.strike_price /1000;
% option prices = (bid + ask)/2
IVSPX.price = (IVSPX.best_bid + IVSPX.best_offer)./2;
    
% VIX data (SECID = 117801)

%% term strucutre of risk-free rate 
df = datastore('rf_term_final.csv');
%preview(df)
rf = readall(df);

%% underlying asset price information
ds = datastore('spx_price_final.csv');
%preview(ds)
ds.SelectedVariableNames={'secid','date','close'};
price_spx = readall(ds);
%%
dvix=datastore('VIX.csv');
dvix.SelectedVariableNames={'date','close'};
VIX_answer=readall(dvix);
%% data construction
% obtain all dates
alldate = unique(IVSPX.date);
%%
% create a space for results [date,oi_var]
% for one-month (21-day) option-implied variance
VIX=zeros(1,1);
SVIX=zeros(1,1);
var=zeros(1,1);
svar=zeros(1,1);
finalsvar=zeros(1,1);
corr=zeros(1,1);
for i = 1:length(alldate)
    % 找到日期 i 的資料 data
    data = IVSPX(IVSPX.date == alldate(i),:);
    % 找到日期 i 的標的資產價格
    St = price_spx.close(price_spx.date == alldate(i));  
    % 找到日期 i 的無風險利率
    %rff = rf(rf.)
    
   %%
    
    % 計算到期天數 maturity
    % business day的計算方式
    data.maturity = days252bus(datetime(data.date,'ConvertFrom','yyyymmdd'),datetime(data.exdate,'ConvertFrom','yyyymmdd'));
    % 刪除到期日 < 5 (business days)
    data = data(data.maturity >= 5 & data.maturity <=60,:);
    allmat = unique(data.maturity); % collect all maturity in each day
    if or(any(allmat(:)==21),~(any(allmat(:)-21<0&any(allmat(:)-21>0))))%如果剛好有21 或是mat太少 比方說只有24 43 51 沒有21以下
        % obtain risk-free rate
        check=0;
        mat = 21;
        rf_d = rf(rf.date == alldate(i),:);
        r = interp1(rf_d.days,rf_d.rate,21,'linear','extrap');
        % obtain maturity
        T = mat/252;

        [~,index] = min(abs(allmat-21));
        data = data(data.maturity == allmat(index),:);
        data = sortrows(data,'strike_price');
        %算F
        strikep=unique(data.strike_price);
        temp_c_p=zeros(1,4);
        for p=1:length(strikep)
            temp_c_p(p,1)=strikep(p);
            temp=data.price(data.strike_price==strikep(p));
            if  isequal(size(temp),[1,1])%有可能只有C或只有P
                temp_c_p(p,4)=99999;%造一個超大值 不要取到那個
                continue
            end
            temp_c_p(p,2)=temp(1);%C
            temp_c_p(p,3)=temp(2);%P
            temp_c_p(p,4)=abs(temp(1)-temp(2));
        end
        [~,index_for_argmin_strikeprice]=min(temp_c_p(:,4));
        ATMstrike=strikep(index_for_argmin_strikeprice);
        ATMForward=ATMstrike+exp(r*T)*(temp_c_p(index_for_argmin_strikeprice,2)-temp_c_p(index_for_argmin_strikeprice,3));
        midquote=data.price(data.strike_price==ATMstrike);
        % 保留OTM call
        % 保留OTM put
        data = data(strcmp(data.cp_flag,'C') & (data.strike_price>St)|strcmp(data.cp_flag,'P') & (data.strike_price<St),:);

        check_for_bid=data.best_bid==0;%把bid=0的地方標記為1
        check_for_continuous=check_for_bid(1:end-1)+check_for_bid(2:end);%連續兩個bid=0 相加會得2 所以找出2所包住的範圍即為所求
        check_index=find(check_for_continuous==2);%把2的index列出 包住的部分即為所求。
        [~,mid]=max(data.volume);%找一個中點 bidprice一定不為0的地方

        temp_max=min(find(check_index>mid));
        temp_min=max(find(check_index<mid));
        if(isempty(temp_max))
            index_max=length(check_for_bid);
        else
            index_max=check_index(temp_max)-1;
        end
        if(isempty(temp_min))
            index_min=1;
        else
            index_min=check_index(temp_min)+2;
        end
        data=data(index_min:index_max,:);
        data=data(data.best_bid~=0,:);

        % add new variable 'deltaK'
        data.deltaK = [0; data.strike_price(2:end) - data.strike_price(1:end-1)];
        %刪除變數Var11
        %data = removevars(data, 'Var11');

        %%
        data_price_adjust=data.price;
        strike_index=find(data.strike_price==ATMstrike);
        data_price_adjust(strike_index)=mean(midquote);
        data.sum_part=(data.deltaK./power(data.strike_price,2))*exp(0.01*r*T).*data_price_adjust; %用.*和./ 而不是直接* / 因為matlab這樣是矩陣乘法
        data.sum_part_SVIX=(data.deltaK.*data_price_adjust)*exp(0.01*r*T);
        var(i,1)=(2/T)*sum(data.sum_part)-(1/T)*power((ATMForward/ATMstrike-1),2);
        svar(i,1)=(2/(T*power(ATMForward,2)))*sum(data.sum_part_SVIX)-(1/T)*power((ATMForward/ATMstrike-1),2);
        VIX(i,1)=100*sqrt(var(i,1)*T*525600/43200);
        SVIX(i,1)=100*sqrt(svar(i,1)*T*525600/43200);
        var(i,1)=var(i,1)*T*525600/43200;
        finalsvar(i,1)=svar(i,1)*T*525600/43200;
    else
        check=1;
        mat = zeros(2,1);
        T=zeros(2,1);
        allmat_1 = allmat(allmat<21);
        allmat_2 = allmat(allmat>21);% 大於才對?
        [~,index] = min(abs(allmat_1-21));
        mat(1) = allmat_1(index);
        [~,index] = min(abs(allmat_2-21));
        mat(2) = allmat_2(index);
        for j = 1:2
            data2 = data(data.maturity == mat(j),:);
            data2 = sortrows(data2,'strike_price');
            
            % obtain risk-free rate
            rf_d = rf(rf.date == alldate(i),:);
            r = interp1(rf_d.days,rf_d.rate,mat(j),'linear','extrap');
            % obtain maturity
            T(j) = mat(j)/252;
            %算F
            strikep=unique(data2.strike_price);
            temp_c_p=zeros(1,4);
            for p=1:length(strikep)
                temp_c_p(p,1)=strikep(p);
                temp=data2.price(data2.strike_price==strikep(p));
                if  isequal(size(temp),[1,1])%有可能只有C或只有P
                    temp_c_p(p,4)=99999;%造一個超大值 不要取到那個
                    continue
                end
                temp_c_p(p,2)=temp(1);%C
                temp_c_p(p,3)=temp(2);%P
                temp_c_p(p,4)=abs(temp(1)-temp(2));
            end
            [~,index_for_argmin_strikeprice]=min(temp_c_p(:,4));
            ATMstrike=strikep(index_for_argmin_strikeprice);
            ATMForward=ATMstrike+exp(r*T(j))*(temp_c_p(index_for_argmin_strikeprice,2)-temp_c_p(index_for_argmin_strikeprice,3));
            midquote=data2.price(data2.strike_price==ATMstrike);

            % 保留OTM call
            % 保留OTM put
            data2 = data2(strcmp(data2.cp_flag,'C') & (data2.strike_price>St)|strcmp(data2.cp_flag,'P') & (data2.strike_price<St),:);
    
            check_for_bid=data2.best_bid==0;%把bid=0的地方標記為1
            check_for_continuous=check_for_bid(1:end-1)+check_for_bid(2:end);%連續兩個bid=0 相加會得2 所以找出2所包住的範圍即為所求
            check_index=find(check_for_continuous==2);%把2的index列出 包住的部分即為所求。
            [~,mid]=max(data2.volume);%找一個中點 bidprice一定不為0的地方 這裡用volume最多地方作為中點 確保2包的位置含有這個中點
    
            temp_max=min(find(check_index>mid));
            temp_min=max(find(check_index<mid));
            if(isempty(temp_max))
                index_max=length(check_for_bid);
            else
                index_max=check_index(temp_max)-1;
            end
            if(isempty(temp_min))
                index_min=1;
            else
                index_min=check_index(temp_min)+2;
            end
            forcheck1=data2;
            data2=data2(index_min:index_max,:);
            forcheck2=data2;
            data2=data2(data2.best_bid~=0,:);
            forcheck3=data2;
            data2.deltaK = [0; data2.strike_price(2:end) - data2.strike_price(1:end-1)]; 
            %%
            data_price_adjust=data2.price;
            strike_index=find(data2.strike_price==ATMstrike);
            data_price_adjust(strike_index)=mean(midquote);
            data2.sum_part=(data2.deltaK./power(data2.strike_price,2))*exp(0.01*r*T(j)).*data_price_adjust; %用.*和./ 而不是直接* / 因為matlab這樣是矩陣乘法
            data2.sum_part_SVIX=(data2.deltaK.*data_price_adjust)*exp(0.01*r*T(j));
            var(i,j)=(2/T(j))*sum(data2.sum_part)-(1/T(j))*power((ATMForward/ATMstrike-1),2);
            svar(i,j)=(2/(T(j)*power(ATMForward,2)))*sum(data2.sum_part_SVIX)-(1/T(j))*power((ATMForward/ATMstrike-1),2);

            % continue...

        end
        wnear=(mat(2)-21)/(mat(2)-mat(1));

        VIX(i,1)=100*sqrt((T(1)*var(i,1)*wnear+T(2)*var(i,2)*(1-wnear))*525600/43200);
        SVIX(i,1)=100*sqrt((T(1)*svar(i,1)*wnear+T(2)*svar(i,2)*(1-wnear))*525600/43200);
        var(i,1)=(T(1)*var(i,1)*wnear+T(2)*var(i,2)*(1-wnear))*525600/43200;
        finalsvar(i,1)=(T(1)*svar(i,1)*wnear+T(2)*svar(i,2)*(1-wnear))*525600/43200;
        %oi_var = interp1(mat,var,'linear','extrap');

    end
    %{
    if(i>1)
        corrmatrix=corrcoef(VIX,VIX_answer.close(1:i));
        corr(i,1)=corrmatrix(1,2);
    else
        corr(i,1)=1;
    end
    %}
end
 




    


