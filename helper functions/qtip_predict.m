function signal = qtip_predict(model,bten)
% function qtip_predict()
%
% returns the signal values as implied by the QTI model given the model
% parameters and the b-tensors
%
if size(model,1) == 28
    model = model';
end

bt2 =  convert_1x6_to_1x21(bten);
signal = model(1) * exp(- bten * model(2:7)' + 0.5 * bt2 * model(8:28)');

end