%searchit1.m
%
%v1.3; 25.6.2021
%
%A toy model looking at "creative" search and control processes in a 2-d 
%item space with the Euclidean distance between items reflecting the
%strength of their association
%
%Lines starting %@@@ are where things that could be usefully added in
%future versions
%
%Free parameters to play with are marked by comment on previous line that
%looks like this %>>>>>>>

%commands which set up the script
clc;
clear variables;
rng('default'); %sets the random number generator to its default starting posotion

%simulation parameters
nsims=200; %number of simulations
nitems=100; %number of items in the search space

%item parameters (the semantic network mapping)
%NB two dimensions only
ndims=2;
%the following parameters set up a mean and a variance-covariance matrix
%for generating the semantic space of items
mu=[0 0]; %this is the mean of item values on the 2 dimensions
itemsd=[1 1]; %[sqrt(0.5) sqrt(0.5)]; %[1 1]; %sd on each dimension
itemvar=itemsd.^2; %variance on each dimension
itemcov=0; %covariance between dimensions 1 and 2 set to zero (an assumption)
sigma=[itemvar(1) itemcov; itemcov itemvar(2)]; %this is the variance covariance matrix

%random walk parameters
mystart=[0 0];  %search starts by deafult at the centre of the space; 
%the value of the next parameter is currently chosen to allow plenty of time to
%generally find a response
nsteps=10000; %1500 for testing the model; %this is the total number of steps for the random walk in each simulation
%>>>>>>>
%if you push the next parameters too large then you will get lots of no
%response trials eg if you use 0.5, 0.5
stepsize=[0.05 0.05]; %the max size of movement for each dimension on a single step
%the next parameter controls the duration of "mind wandering"
%>>>>>>>
walkfor=50; %20 50 100. The number of steps of random walk before trying to retrieve an item

%Noisy choice parameter for retrieval of items
%tau is the so-called temperature parameter for noisy choice, controls noisiness
%as values go increasingly below 1 the choice becomes equally likely even if the 
%basis (distance of items from current position) for the choice is quite varied
%higher tau values increase the focus of the choice on the closest items to
%a greater extent
%>>>>>>>
tau = 2; %0.2 0.5 1 2 5
%>>>>>>>
closeto=1.0; %0.5 0.75 1. This sets which items can enter choice of retrieval (i.e. any closer than this Euc dist)
smallval=0.01; %this just prevents overflow errors when using the softmax choice rule

%parameters affecting decision whether a retrieved item should be given as a response
respmethod=1; %use the threshold on Euc distance from centre
%respmethod=2; %this is a null choice that is used for simulating a fluency task
%if you change the next parameter it trivially has a direct effect on the
%response choice and time to reach it
respthreshbase=2; %this sets the distance from the centroid which a retrieved item must be in order to be given as a response
%>>>>>>>
threshdrop= -0.01; %0 -0.01 -0.02 can be used here; this sets how much respthresh decreases per number of timesteps set by threstime
threshtime=100;

%initialise key arrays for recording results of simulations
%they are initialised with not a number, NaN
%so we can filter cases where the simulation does not produce a value
myresp=NaN(nsims,1); %this will store the response number of each item in each sim
myresptime=NaN(nsims,1); %this will store the response time in each sim
myresppos=NaN(nsims,ndims); %the position of the response items in the 2-d space
myretrieves=NaN(nsims, nitems); %records the items retrieved
myretrievaltimes =NaN(nsims, nitems); %records the times when the items were retrieved
edistfromcent=NaN(nsims,1); %this will store the Eucidean dist from centre of the responses
num_retrieved=NaN(nsims,1); %this will store the number of retrievals made in each simulation
unique_retrieves=NaN(nsims,1); %this will store the number of unique retrievals made in each simulation

%plotting and printing controls
%simbysim=1 gives a summary print to screen at end of each sim; 0 otherwise
simbysim=-1;
disp('Results will display at the end of the simulations.');
while simbysim~=0 && simbysim~=1
    simbysim=input('Do you ALSO want a summary result after each simulation (1=yes; 0=no) + <Enter> ');
end
plot2show=1;        %which simulation to plot for (>nsims=no plot)
print2show=20;      %controls printing of onscreen messages, ie every print2show sims
disp('Simulating...');
%we imagine each sim as a unique participant
for s=1:nsims
    
    %display sim counter
    if mod(s, print2show)==0 && simbysim==0
        disp(['Completed sim number ' num2str(s)]);
    end
    
    %create the values on two dimensions for the stimulus exemplars
    %this will be randomly varying for each simulation but the same for
    %each run of the code. The starting points are set far apart by using the 10000 multiplier
    randseed=s*10000;
    rng(randseed,'twister'); %this ensures that the random items in a specific simulation are the same each time code is run
    itemvals=mvnrnd(mu, sigma, nitems); %using a multivariate normal random distribution
    
    %plot out an example item graph
    %for a particular simulation (when sim number, s=plot2show)
    if s==plot2show
        figure;
        plot(itemvals(:,1),itemvals(:,2),'ok');
        title(['Position of items in simulation #' num2str(s)]);
        xlabel('Dimension 1');
        ylabel('Dimension 2');
    end
    
    %run a random walk for a fixed amount of time 
    currpos=mystart; 
    %compute a random movement -1 0 or 1 for each dimension separately
    %so 9 possible directions for moves, all equiprobable 0 0, 0 +1, 0 -1 etc
    %unique random set of moves for each simulation
    mymoves=randi([-1 1],nsteps,ndims);
    endkflag=0; %this is just a flag to control when the walk needs to finish early
    numretrieved=0;
    respthresh=respthreshbase;
    for k=1:nsteps
        
        if mod(k,threshtime)==0
            respthresh=respthresh+threshdrop; %decreases threshold because threshdrop is -ve
        end
        
        %do the walk
        %NB multiply the move diections in each direction by the step size
        %in that direction
        currpos=currpos+stepsize.*(mymoves(k,:));
        
        %now pause the walk to enact retrieval
        if mod(k,walkfor)==0
            %compute euclidean distances from current position to all the items
            eucdist=sqrt((currpos(1)-itemvals(:,1)).*(currpos(1)-itemvals(:,1)) + (currpos(2)-itemvals(:,2)).*(currpos(2)-itemvals(:,2)));
             
            %here we can add noisy choice for which item gets retrieved
            %we filter only the relatively close exemplars
            %ie those closer than closeto
            choicefilt=eucdist<closeto;
            if sum(choicefilt)>0
                %this if is needed as we can retrieve an item only if there 
                %is at least 1 within the required range
                numretrieved=numretrieved+1; %counts the number retrieved during this simulation
                
                %NB we use the inverse of the eucdist (1/ED) as the data in
                %the choice function, but very small values of ED lead to
                %very large (1/ED) values. This can cause Inf values for
                %probchoice. So all the "close enough" eucdist values are set
                %to be >= smallval (0.01) to prevent this happening
                eucdist(choicefilt)=max(eucdist(choicefilt),smallval);
                denom=sum(exp(tau./eucdist(choicefilt))); %this is the denominator of the choice function
                probchoice=zeros(nitems,1); %set all probs to zero first
                probchoice(choicefilt)=exp(tau./eucdist(choicefilt))./denom; %choice probs just for the close items
                
                %next we have to execute the choice randomly according
                %the values computed above in probchoice
                %this uses cumulative sums and the uniform random number
                %generator, rand
                cumprobs=cumsum(probchoice(choicefilt));
                mychoice=find(rand<cumprobs,1); %find the first value in the cumulative probs greater than the rand num
                retrieveditem=max(find(choicefilt,mychoice)); %find the item number in the filtered item list for mychoice
                
                %we store the retrieved items here and the retrieval times
                %per simulation
                myretrieves(s,numretrieved)=retrieveditem;
                myretrievaltimes(s,numretrieved)=k;
                edistfromcent(s)=sqrt((itemvals(retrieveditem,1)-mu(1)).*(itemvals(retrieveditem,1)-mu(1))+ ...
                    (itemvals(retrieveditem,2)-mu(2)).*(itemvals(retrieveditem,2)-mu(2))); %this is the distance from the centroid
                %and here we move the current position to the position of
                %the retrieved item
                currpos=itemvals(retrieveditem,:);
                %here we decide if ithe item is obscure enough to be our response
                if respmethod==1
                    %the simple decision rule here is that we have a threshold distance from
                    %the centroid (0,0) of the item space and if that threshold
                    %is exceeded by the retrieved item, then we select that item
                    if edistfromcent(s) > respthresh
                        %the next lines of code are useful for
                        %visualising the effect of tau in a simulation
                        if s==plot2show
                            figure;
                            plot(eucdist(choicefilt), probchoice(choicefilt), 'xk');
                            xlabel('Euclidean distance from current position');
                            ylabel('Probability of retrieval');
                            title({'Retrievable items just before selection made', ['Simulation #' num2str(s) ' tau= ' num2str(tau)]});
                        end
                        endkflag=1; 
                        %records response selection
                        myresp(s)=retrieveditem;
                        myresppos(s,:)=itemvals(retrieveditem,:);
                        myresptime(s,1)=k;
                    else
                        %item not chosen as response
                        endkflag=0;
                    end
                elseif respmethod==2
                    %currently this turns off the response decision
                    if 0 %note this condition will never be satisfied
                        endkflag=1;  %#ok<UNRCH>
                        %record response selection for code completeness
                        myresp(s)=retrieveditem;
                        myresppos(s,:)=itemvals(retrieveditem,:);
                        myresptime(s,1)=k;
                    else
                        %item not chosen as response
                        endkflag=0;
                    end
                    %@@@could add other decision rule (respmethod) choices here
                    %end of respmethod if loop
                end
                
            else
                %do nothing if none are close enough for retrieval choice
            end
            %end of if loop controlling whether or not any items are close enough
            
        %end of if loop controlling when we pause in walk for retrieval
        end
        
        if endkflag==1
            break;  %to ensure no more walking on this sim
        end
            
        
        %end of random walk steps, k loop
    end
    
    %@@@in a full model we would need to choose a response even if none of the retrieved items have passed our
    %response decision criterion. 

    %compute number of items retrieved in this sim
    num_retrieved(s)=sum(~isnan(myretrieves(s,:)));
    unique_retrieves(s)=size(unique(myretrieves(s,~isnan(myretrieves(s,:)))),2);
    
    %show some information simulation by simultion if requested
    if simbysim==1
        disp(['Simulation #' num2str(s)]);
        disp(['Number of retrievals = ' num2str(num_retrieved(s))]);
        disp(['Number of unique items retrieved = ' num2str(unique_retrieves(s))]);
        disp('Retrieved items and their retrieval times');
        showit=[myretrieves(s,~isnan(myretrieves(s,:)))', myretrievaltimes(s,~isnan(myretrievaltimes(s,:)))'];
        disp(showit);
        disp(['Response made= ' num2str(myresp(s))]);
        disp(['Position of resp made= ' num2str(myresppos(s,:))]);
        disp(['Euc. dist. of resp made (or last item retrieved) from centre = ' num2str(edistfromcent(s))]);
        disp('Hit any key for next simulation')
        disp(' ');
        pause;
    end
            
%end of simulation, s loop
end

%now the results across all simulations
%@@@could display key parameter values that are being used here

%draw where the responses lie in the 2-d space
figure;
plot(myresppos(:,1),myresppos(:,2),'ok');
xlabel('Dimension 1');
ylabel('Dimension 2');
title('Positions of items given as responses across simulations');

%draw histogram of response times
figure;
histogram(myresptime);
xlabel('Response times (in walk steps)');
ylabel('Frequency');
title('Response time distribution');
%@@@could add other histograms

%key results as text displayed on screen
disp(' ');
disp('Retrieval info across all sims');
disp('------------------------------');
disp(['The number of simulations where at least 1 item was retrieved = ' num2str(sum(~isnan(unique_retrieves))) ' / ' num2str(nsims)]);
disp(['The mean number of retrieved items    = ' num2str(mean(num_retrieved,'omitnan'))]);
disp(['The s.d. of number of retrieved items = ' num2str(std(num_retrieved,'omitnan'))]);
disp(['The mean number of unique retrieved items    = ' num2str(mean(unique_retrieves,'omitnan'))]);
disp(['The s.d. of number of unique retrieved items = ' num2str(std(unique_retrieves,'omitnan'))]);


disp(' ');
disp('Response info across all sims');
disp('-----------------------------');
disp(['The number of simulations where a response was made = ' num2str(sum(~isnan(myresp))) ' / ' num2str(nsims)]);
disp(['The mean position of the item responses        = ' num2str(mean(myresppos,'omitnan'))]);
disp(['The s.d. of the position of the item responses = ' num2str(std(myresppos,'omitnan'))]);
respdist(:,1)=sqrt((myresppos(:,1)-mu(1)).*(myresppos(:,1)-mu(1)) + (myresppos(:,2)-mu(2)).*(myresppos(:,2)-mu(2)));
disp(['The mean Euclidean distance of the item responses from the centroid = ' num2str(mean(respdist,'omitnan'))]);
disp(['The s.d. Euclidean distance of the item responses from the centroid = ' num2str(std(respdist,'omitnan'))]);
disp(['The mean times of item responses (in random walk steps)        = ' num2str(mean(myresptime,'omitnan'))]);
disp(['The s.d. of the times of item responses (in random walk steps) = ' num2str(std(myresptime,'omitnan'))]);
