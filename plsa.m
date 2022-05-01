function [prob_term_topic, prob_topic_doc, lls,time,errorW,errorH,Terror,iter] = plsa(termDocMatrix,numTopic,iter)
tstart=tic;
time=zeros(10000,1);
errorW=time;
errorH=time;
Terror=time;
% record error
[numTerm, numDoc] = size(termDocMatrix);



% Fit a plsa model from a given term-document matrix



prob_term_topic = rand(numTerm, numTopic); % p(term | topic)
for i = 1:numTopic
	prob_term_topic(:, i) = prob_term_topic(:, i) / sum(prob_term_topic(:, i));
end

prob_topic_doc = rand(numTopic, numDoc);   % p(topic | doc)
for i = 1:numDoc
	prob_topic_doc(:, i) = prob_topic_doc(:, i) / sum(prob_topic_doc(:, i));
end
% for l=1:numTopic
% TW=prob_term_topic.^(1/2)-trW.^(1/2);
% errorW(1)=errorW(1)+(sum((TW(:,l)).^2)/2)^(1/2)/numTopic;
% end
% for l=1:numDoc
% TH=prob_topic_doc.^(1/2)-trH.^(1/2);
% errorH(1)=errorH(1)+(sum((TH(:,l)).^2)/2)^(1/2)/numDoc;
% end
V0=prob_term_topic*prob_topic_doc;
TV=V0.^(1/2)-termDocMatrix.^(1/2);
for l=1:numDoc
Terror(1)=Terror(1)+(sum((TV(:,l)).^2)/2)^(1/2)/numDoc;
end

prob_term_doc = zeros(numTerm, numDoc);
for d = 1:numDoc
	for z = 1:numTopic
		prob_term_doc(:, d) = prob_term_doc(:, d) + ...
		prob_topic_doc(z, d) .* prob_term_topic(:, z);
	end
	assert(sum(prob_term_doc(:, d)) - 1.0 < 1e-6);
end

prob_topic_term_doc = cell(numTopic, 1);   % p(topic | doc, term)
for z = 1 : numTopic
	prob_topic_term_doc{z} = zeros(numTerm, numDoc);
end


lls = []; % maximum log-likelihood estimations

for i = 1 : iter
	disp('E-step...');
	for d = 1:numDoc
		%fprintf('processing doc %d\n', d);
		w = find(termDocMatrix(:, d));
		for z = 1:numTopic
			prob_topic_term_doc{z}(w, d) = prob_topic_doc(z, d) .* prob_term_topic(w, z) ./ prob_term_doc(w, d);
		end
	end
	
	disp('M-step..');
	disp('Update p(topic | doc)...');
	for d = 1:numDoc
		w = find(termDocMatrix(:, d));
		for z = 1:numTopic
% 			prob_topic_doc(z, d) = max(sum(termDocMatrix(w, d) .* prob_topic_term_doc{z}(w, d))-0.1/numTopic,0);
            prob_topic_doc(z, d) = max(sum(termDocMatrix(w, d) .* prob_topic_term_doc{z}(w, d)),0);
		end
		prob_topic_doc(:, d) = prob_topic_doc(:, d) / sum(prob_topic_doc(:, d));
	end
	disp('Update p(word | topic)...');
	for z = 1:numTopic
		for w = 1:numTerm
			d = find(termDocMatrix(w, :));
% 			prob_word_topic(w, z) = max(sum(termDocMatrix(w, d) .* prob_topic_term_doc{z}(w, d))-0.1/numTerm,0);
            prob_word_topic(w, z) = max(sum(termDocMatrix(w, d) .* prob_topic_term_doc{z}(w, d)),0);
		end
		prob_word_topic(:, z) = prob_word_topic(:,z) / sum(prob_word_topic(:,z));
	end
	
	% calculate likelihood and update p(term, doc)
	fprintf('Iteration %d\n', i);
	disp('Calculate maximum likelihood...');
	ll = 0;
	for d = 1:numDoc
		prob_term_doc(:, d) = 0;
		for z = 1:numTopic
			prob_term_doc(:, d) = prob_term_doc(:, d) + ...
			prob_topic_doc(z, d) .* prob_term_topic(:, z);
		end
		assert(sum(prob_term_doc(:, d)) - 1.0 < 1e-6);
		w = find(termDocMatrix(:, d));
		ll = ll + sum(termDocMatrix(w, d) .* log(prob_term_doc(w, d)));
	end
	fprintf('likelihood: %f\n', ll);
	lls= [lls;ll];
    
% for l=1:numTopic
% TW=prob_term_topic.^(1/2)-trW.^(1/2);
% errorW(i+1)=errorW(i+1)+(sum((TW(:,l)).^2)/2)^(1/2)/numTopic;
% end
% 
% for l=1:numDoc
% TH=prob_topic_doc.^(1/2)-trH.^(1/2);
% errorH(i+1)=errorH(i+1)+(sum((TH(:,l)).^2)/2)^(1/2)/numDoc;
% end
V0=prob_term_topic*prob_topic_doc;
TV=V0.^(1/2)-termDocMatrix.^(1/2);
for l=1:numDoc
Terror(i+1)=Terror(i+1)+(sum((TV(:,l)).^2)/2)^(1/2)/numDoc;
end

tend=toc(tstart);
time(i+1)=tend;
    
    
end
% save model.mat prob_topic_doc prob_word_topic
end
