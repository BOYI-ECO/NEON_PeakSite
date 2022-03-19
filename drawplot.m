function []=drawplot(xx,xreal,range,mm,nn,yla)
% xx: parameter samples obtained with inverse methods
% xreal: real parameter values
% range: prior ranges for parameters
% mm,nn: subplot(mm,nn,i)
% yla: the name of each parameters
titlestr={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)'...
    ,'(o)','(p)','(q)','(r)','(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};

mn = mm*nn;
[a,b]=size(xx);
if a>b
    xx=xx';
end
[a,~]=size(xx);
if nargin<6
    for i=1:a
        yla{i}=['m',num2str(i)];
    end
end

mn = mm*nn;
if a <= mn
    figure('Color',[1 1 1]);
    for i=1:a
        subplot(mm,nn,i,'FontWeight','bold','FontSize',12)
        plot(xx(i,:),'MarkerSize',15,'Marker','.','Linestyle','None');hold on;
        if isempty(xreal) ==0
            plot([0 size(xx,2)],[xreal(i) xreal(i)],'LineWidth',2,'Color',[1 0 0]);
        end
        hold off;
        axis([0 size(xx,2) range(i,1) range(i,2)]);
        xlabel('Number of model evaluations','FontWeight','bold','FontSize',10);
        ylabel(yla{i},'FontWeight','bold','FontSize',10)
        %         text(size(xx,2)*0.85,range(i,2)-(range(i,2)-range(i,1))*0.18,titlestr{i},'EdgeColor',[1 1 1],...
        %             'BackgroundColor',[1 1 1],'FontWeight','bold','FontSize',12);
    end
    legend('Parameter samples','True value');
end

if a > mn
    e = ceil(a/mn)-1;
    f = mod(a,mn);
    if f==0
        ee=e+1;
    else
        ee=e;
    end
    for j=1:ee
        figure('Color',[1 1 1]);
        for i=1:mn
            subplot(mm,nn,i,'FontWeight','bold','FontSize',12)
            plot(xx(j*mn-mn+i,:),'MarkerSize',15,'Marker','.','Linestyle','None');hold on;
            if isempty(xreal)==0
                plot([0 size(xx,2)],[xreal(j*mn-mn+i) xreal(j*mn-mn+i)],'LineWidth',2,'Color',[1 0 0]);
            end;
            hold off;
            axis([0 size(xx,2) range(j*mn-mn+i,1) range(j*mn-mn+i,2)]);
            xlabel('Number of model evaluations','FontWeight','bold','FontSize',10);
            ylabel(yla{j*mn-mn+i},'FontWeight','bold','FontSize',10)
            %             text(size(xx,2)*0.85,range(j*mn-mn+i,2)-(range(j*mn-mn+i,2)-range(j*mn-mn+i,1))*0.18,titlestr{i},'EdgeColor',[1 1 1],...
            %                 'BackgroundColor',[1 1 1],'FontWeight','bold','FontSize',12);
        end
    end
    if f>0
        figure('Color',[1 1 1]);
        for i=1:f
            subplot(mm,nn,i,'FontWeight','bold','FontSize',12)
            plot(xx(e*mn+i,:),'MarkerSize',15,'Marker','.','Linestyle','None');hold on;
            if isempty(xreal)==0
                plot([0 size(xx,2)],[xreal(e*mn+i) xreal(e*mn+i)],'LineWidth',2,'Color',[1 0 0]);
            end
            hold off;
            axis([0 size(xx,2) range(e*mn+i,1) range(e*mn+i,2)]);
            xlabel('Number of model evaluations','FontWeight','bold','FontSize',10);
            ylabel(yla{i+e*mn},'FontWeight','bold','FontSize',10)
            %             text(size(xx,2)*0.85,range(e*mn+i,2)-(range(e*mn+i,2)-range(e*mn+i,1))*0.18,titlestr{i},'EdgeColor',[1 1 1],...
            %                 'BackgroundColor',[1 1 1],'FontWeight','bold','FontSize',12);
        end
    end
end
if isempty(xreal)==0
    legend('Parameter samples','True value');
else
    legend('Parameter samples');
end
end


