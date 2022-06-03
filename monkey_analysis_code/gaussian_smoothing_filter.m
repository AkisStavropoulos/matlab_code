%% filtering data with gaussian filter

%% my signal
Nt = 1000;
data = normrnd(0,1,[1 Nt]);

%% my filter
L = 30; sig = 10;       % L is not playing an important role if L >> sigma. L has to be big enough (number of samples taken)
t = linspace(-L,L,2*L+1);

my_filter = exp(-(t - 0).^2/(2*sig^2));
my_filter = my_filter/sum(my_filter);

%% my filtered signal
my_filtered_data = conv(data,my_filter,'same');

%% plot
figure;
subplot(2,2,1);
plot(1:Nt,data);
title('given fuction');

subplot(2,2,2);
plot(t,my_filter);
title('filter applied'); xlabel('L = 95% of distribution, \sigma = 60%');

subplot(2,2,3);
plot(1:Nt,my_filtered_data);
title('filtered data');