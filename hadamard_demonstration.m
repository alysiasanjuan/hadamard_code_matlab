clc;
while 1 
    disp('---------------------------------------------------'); disp(' ');
    % Define N, K, and the message
    N = input('Enter the order of Hadamard matrix (N): '); % codeword bits per block
    if rem(log2(N),1) ~= 0 || N == 1 || N == 2
        disp(' '); disp(['⛔ Invalid input, N = ', num2str(N), ' is not a power of 4 ⛔']); 
        continue;
    end
    K = nextpow2(N); % message bits per block1
    m = randi([0 1], 1, K);%[1 1 1]; %
    m0 = [m zeros(1, N-K)];
    err = 2^(K-2)-1;
    disp(' ');
    disp(['Demonstration of Hadamard(', num2str(N), ', ', num2str(K), ')']); 
    disp(['Can correct up to ', num2str(err), ' errors']); 
    err = input('Enter the number of errors to be introduced: ');
    disp(' ');
    
    % Define the Hadamard matrix and Generator matrix
    H = hadamard(N);
    Hx = H;
    Hx(Hx==1)=0;
    Hx(Hx==-1)=1;
    %Hx(Hx==-1)=0;
    G = dec2bin(0:N-1,K) - 48;
    Gp = G';
    
    
    % Encode message
    code = mod(m * Gp, 2);
    disp(['Original message: ', num2str(m)]); 
    disp([ 'Encoded to:   ', num2str(code)]); disp(' ');
    
    % Introduce noise
    noisy = code;
    for i = 1:err
        errindex = randi([1 N], 1, 1);
        noisy(errindex) = mod(noisy(errindex) + 1, 2); % randomize noise
    end
    %noisy(3) = mod(noisy(3) + 1, 2); % make noise random lol
    disp(['Received message: ', num2str(noisy)]);
    
    % Decode message
    syndrome = mod(noisy * Hx, 2);
    if sum(syndrome) == 0
        disp('✅ No error detected ✅');
        [tf, index]=ismember(noisy,Hx,'rows');
        decode = dec2bin(index-1) - 48;
        disp(' '); disp(['Decoded to:   ', num2str(decode)]); disp(' ');
    else
        x = abs((2*noisy - H(1,:))*H);
        [M,I] = max(x,[],2,'linear'); %grabs the index of max as that is the closest to codewords
        U=unique(x);
        u = U(1<histc(x,unique(x)));
        if any(u(:) == M)
	        disp("⛔ Too many errors, impossible to decode ⛔");
        else
            fixed = mod(2*Hx(1,:)-Hx(I,:),2); %bro it looks so complicated but i was jus tryna transform -1 nd 1 to 0 nd 1 ugh
            disp('⚠️ Received message not a codeword️ ⚠️');
            disp(['Closest codeword: ', num2str(fixed)]);
            [tf, index]=ismember(fixed,Hx,'rows');
            decode = dec2bin(index-1) - 48;
            disp(' '); disp(['Decoded to:   ', num2str(decode)]); disp(' ');
        end
    end
    break;
end    

