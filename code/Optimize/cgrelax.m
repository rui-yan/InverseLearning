function [f,x] = cgrelax(func,n,acc,maxfn,dfpred,x0,param,f_is_constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stand alone conjugate gradient relaxiation program
%%
%%    Based on IMSL's U2CGG.F
%%
%%  rewritten in C by Dongyi Liao at MIT 1999
%%  rewritten in matlab by Wei Cai at LLNL 2004
%%
% @ param: func, target function to minimize
% @ param: n, number of variables
% @ param: acc, accuracy, namely the squared sum of residual gradient 
% @ param: maxfn, max # iterations
% @ param: dfpred, rough estimate of the expected reduction in the
% function. Usually used for determining the size of initial change to x
% @ param: x0, initial value
% @ param: param, param for function f
% @ global: print_level, if == 1, print out computation details
% @ returns: f, final function value
% @ returns: x, vector of length n containing the computed solution.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global print_level
global debug

%define
MAXLIN = 5; %% Maximum iterations within a line search
MXFCON = 2; %% Maximum iterations before F can be decreased

%initialize memory
x=zeros(n,1);
g=zeros(n,1);
        
s=zeros(n,1);
rss=zeros(n,1);
rsg=zeros(n,1);
ginit=zeros(n,1);
xopt=zeros(n,1);
gopt=zeros(n,1);

    
ddspln=0;
gamden=0;

iterc=0;
iterfm=0;
iterrs=0;

x=x0;
%[f,g]=feval(func,n,x);
[f,g]=cg_fchem(x,param,f_is_constant);

ncalls=1;
if strcmp(debug,'on')
    disp('relax: 1st potential call finished.');
end

s=-g;
gsum=sum(g.*g);
    
if((gsum)<=acc)     % change from gsum<=acc
    if strcmp(debug,'on')
        disp('converged');
    end
    return;
end

gnew=-gsum;
fmin=f;
gsqrd=gsum;
nfopt=ncalls;

xopt=x;
gopt=g;

dfpr=dfpred;
stmin=dfpred/gsqrd;

reinsert=0;

while(1) 
%% BEGIN THE ITERATION

    if(~reinsert) %% A MECHANISM TO IMPLEMENT RETRY
        
    iterc=iterc+1;
    %% STORE THE INITIAL FUNCTION VALUE AND GRADIENT, CALCULATE THE
    %% INITIAL DIRECTIONAL DERIVATIVE, AND BRANCH IF ITS VALUE IS NOT
    %% NEGATIVE. SET SBOUND TO MINUS ONE TO INDICATE THAT A BOUND ON THE
    %% STEP IS NOT KNOWN YET, AND SET NFBEG TO THE CURRENT VALUE OF
    %% NCALLS. THE PARAMETER IRETRY SHOWS THE NUMBER OF ATTEMPTS AT
    %% SATISFYING THE BETA CONDITION.
    finit=f;
    ginit=g;

    gfirst=sum(s.*g);
    if(gfirst >=0)
        %% SET IER TO INDICATE THAT THE SEARCH DIRECTION IS UPHILL.
        disp('ERROR [cgrelax]: Search direction is uphill.');
        return;
    end
    gmin=gfirst;
    sbound=-1;
    nfbeg=ncalls;
    iretry=-1;
    %% SET STEPCH SO THAT THE INITIAL STEP-LENGTH IS CONSISTENT WITH THE
    %% PREDICTED REDUCTION IN F, SUBJECT TO THE CONDITION THAT IT DOES
    %% NOT EXCEED THE STEP-LENGTH OF THE PREVIOUS ITERATION. LET STMIN
    %% BE THE STEP TO THE LEAST CALCULATED VALUE OF F.
    stepch=min(stmin, abs(dfpr/gfirst));
    stmin=0;
    %% CALL SUBROUTINE FUNCT AT THE VALUE OF X THAT IS DEFINED BY THE
    %% NEW CHANGE TO THE STEP-LENGTH, AND LET THE NEW STEP-LENGTH BE
    %% STEP. THE VARIABLE WORK IS USED AS WORK SPACE.
    
    end %% A MECHANISM TO IMPLEMENT RETRY
    
    while(1)
        if(~reinsert) %% A MECHANISM TO IMPLEMENT RETRY
        
        step=stmin+stepch;
        x=xopt+stepch*s';
        
        temp=max(abs(x-xopt));
        if(temp <= 0)
            %% TERMINATE THE LINE SEARCH IF STEPCH IS EFFECTIVELY ZERO.
            if(ncalls > nfbeg+1 | abs(gmin/gfirst) > 0.2)
                disp('ERROR [cgrelax]: Line search aborted, possible error in gradient.');
                return;
            else
                break;
            end
        end
        ncalls=ncalls+1;
            
        %[f,g]=feval(func,n,x);
        [f,g]=cg_fchem(x,param,f_is_constant);

        %% SET SUM TO G SQUARED. GMIN AND GNEW ARE THE OLD AND THE NEW
        %% DIRECTIONAL DERIVATIVES ALONG THE CURRENT SEARCH DIRECTION.
        gnew=sum(s.*g);
        gsum=sum(g.*g);
            
        %% STORE THE VALUES OF X, F AND G, IF THEY ARE THE BEST THAT
        %% HAVE BEEN CALCULATED SO FAR, AND NOTE G SQUARED AND THE VALUE
        %% OF NCALLS. TEST FOR CONVERGENCE.
        if((f < fmin || (f==fmin && gnew/gmin>=-1)))
            fmin=f;
            gsqrd=gsum;
            nfopt=ncalls;
            xopt=x;
            gopt=g;
        end

        if (print_level==1)
            disp(sprintf('%10d %24.16e %24.16e',iterc,f,gsum));
        end
        
        if(f<=fmin && gsum <= acc)
            if strcmp(debug,'on')
                disp('converged');
            end
            return;
        end

        if(ncalls==maxfn)
            disp('ERROR [cgrelax]: Too many iterations, stop...');
            return;
        end

        %% LET SPLN BE THE QUADRATIC SPLINE THAT INTERPOLATES THE
        %% CALCULATED FUNCTION VALUES AND DIRECTIONAL DERIVATIVES AT THE
        %% POINTS STMIN AND STEP OF THE LINE SEARCH, WHERE THE KNOT OF
        %% THE SPLINE IS AT 0.5*(STMIN+STEP). REVISE STMIN, GMIN AND
        %% SBOUND, AND SET DDSPLN TO THE SECOND DERIVATIVE OF SPLN AT
        %% THE NEW STMIN. HOWEVER, IF FCH IS ZERO, IT IS ASSUMED THAT
        %% THE MAXIMUM ACCURACY IS ALMOST ACHIEVED, SO DDSPLN IS
        %% CALCULATED USING ONLY THE CHANGE IN THE GRADIENT.
        temp=2*(f-fmin)/stepch-gnew-gmin;
        ddspln=(gnew-gmin)/stepch;
        if(ncalls > nfopt) 
            sbound=step;
        else
            if(gmin*gnew <= 0) 
                sbound=stmin;
            end
            stmin=step;
            gmin=gnew;
            stepch=-stepch;
        end
        if(f~=fmin) 
            ddspln=ddspln+(temp+temp)/stepch;
        end
        %% TEST FOR CONVERGENCE OF THE LINE SEARCH, BUT FORCE AT LEAST
        %% TWO STEPS TO BE TAKEN IN ORDER NOT TO LOSE QUADRATIC
        %% TERMINATION.
        if(gmin==0) 
            break;
        end
        if(ncalls >= nfbeg+1)
            if(abs(gmin/gfirst) <=0.2) 
                break;
            end
            %% APPLY THE TEST THAT DEPENDS ON THE PARAMETER MAXLIN.
            %% l_retry:
            if(ncalls >= nfopt+MAXLIN)
                disp('ERROR [cgrelax]: Line search aborted, possible error in gradient.');
                return;
            end
        end
        
        end %% A MECHANISM TO IMPLEMENT RETRY
        
        if(reinsert)%% A MECHANISM TO IMPLEMENT RETRY
            reinsert=0;
            if(ncalls >= nfopt+MAXLIN)
                disp('ERROR [cgrelax]: Line search aborted, possible error in gradient.');
                return;
            end
        end %% A MECHANISM TO IMPLEMENT RETRY
        
        %% SET STEPCH TO THE GREATEST CHANGE TO THE CURRENT VALUE OF STMIN
        %% THAT IS ALLOWED BY THE BOUND ON THE LINE SEARCH. SET GSPLN TO THE
        %% GRADIENT OF THE QUADRATIC SPLINE AT (STMIN+STEPCH). HENCE
        %% CALCULATE THE VALUE OF STEPCH THAT MINIMIZES THE SPLINE FUNCTION,
        %% AND THEN OBTAIN THE NEW FUNCTION AND GRADIENT VECTOR, FOR THIS
        %% VALUE OF THE CHANGE TO THE STEP-LENGTH.
        stepch=0.5*(sbound-stmin);
        if(sbound < -0.5) 
            stepch=9*stmin;
        end
        gspln=gmin+stepch*ddspln;
        if(gmin*gspln<0) 
            stepch=stepch*gmin/(gmin-gspln);
        end
    end %% END OF INNER WHILE LOOP
    
    %% ENSURE THAT F, X AND G ARE OPTIMAL.
    if(ncalls~=nfopt)
        f=fmin;
        x=xopt;
        g=gopt;
    end
    %% CALCULATE THE VALUE OF BETA THAT OCCURS IN THE NEW SEARCH
    %% DIRECTION.
    gsum=sum(g.*ginit);
    beta=(gsqrd-gsum)/(gmin-gfirst);
    %% TEST THAT THE NEW SEARCH DIRECTION CAN BE MADE DOWNHILL. IF IT
    %% CANNOT, THEN MAKE ONE ATTEMPT TO IMPROVE THE ACCURACY OF THE LINE
    %% SEARCH.
    if(abs(beta*gmin) > 0.2*gsqrd)
        iretry=iretry+1;
        if(iretry<=0)
            %% goto l_retry;
            %% disp('need to retry: not implemented yet');
            reinsert = 1;
            continue;
        end
    end
    %% APPLY THE TEST THAT DEPENDS ON THE PARAMETER MXFCON.
    %% SET DFPR TO THE PREDICTED REDUCTION IN F ON THE NEXT ITERATION.
    if(f<finit) 
        iterfm=iterc;
    elseif(iterc >= iterfm+MXFCON)
        disp('ERROR [cgrelax]: Cannot reduce value of F, aborting...');
        return;
    end
    dfpr=stmin*gfirst;
        
    if(iretry>0) %% Restart since we need to retry
        s=-g;
        iterrs=0;
        continue;
    end

    if(iterrs~=0 && (iterc-iterrs)<n && abs(gsum)<0.2*gsqrd)
        %% CALCULATE THE VALUE OF GAMA THAT OCCURS IN THE NEW SEARCH
        %% DIRECTION, AND SET SUM TO A SCALAR PRODUCT FOR THE TEST BELOW.
        %% THE VALUE OF GAMDEN IS SET BY THE RESTART PROCEDURE.
        %% tmp1=tmp2=tmp3=0;
        gama=sum(g.*rsg);
        gsum=sum(g.*rss);
        gama=gama/gamden;
        %% RESTART IF THE NEW SEARCH DIRECTION IS NOT SUFFICIENTLY
        %% DOWNHILL.
        if(abs(beta*gmin+gama*gsum) < 0.2*gsqrd)
            %% CALCULATE THE NEW SEARCH DIRECTION.
            s=-g+beta*s+gama*rss;
            continue;
        end
    end
    %% APPLY THE RESTART PROCEDURE.
    gamden=gmin-gfirst;

    rss=s;
    rsg=g-ginit;
    s=-g+beta*s;
    iterrs=iterc;

end %% END OF WHILE LOOP

%% END OF CGRELAX
