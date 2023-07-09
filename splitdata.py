from re import L


data=[]
count=0
i=0
with open('data_20mb.txt','r') as f:
    for line in f.readlines():
        data.append(line)
        count+=1
    f.close()
data_peruser=int(count/10)
while i <9:
    start=i*data_peruser
    filename='user_'+str(i)+'_20mb.txt'
    with open(filename,'w') as f:
        for j in range(data_peruser):
            f.write(data[start+j])
        
        f.close()
    i+=1
filename='user_'+str(i)+'_20mb.txt'
final_userdata=data_peruser+count%10
start=i*data_peruser
with open(filename,'w') as f:
    for j in range(final_userdata):
        f.write(data[start+j])
