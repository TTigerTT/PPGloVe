vocab_dic={}
for i in range(10):
    file_name='user'+str(i)+'_vocab_20mb.txt'
    with open (file_name ,'r') as f:
        for line in f.readlines():
            line=line.strip()
            line=line.split(' ')
            word=line[0]
            count=int(line[1])

            if word not in vocab_dic:
                vocab_dic[word]=0
            vocab_dic[word]+=count
        f.close()
res=sorted(vocab_dic.items(),key= lambda item:item[1],reverse= True)
print(len(vocab_dic))
with open('mergedvocab_20mb.txt','w') as f:
    for i in res:
        f.write(i[0]+' '+ str(i[1])+'\n')
        
    
    f.close()


