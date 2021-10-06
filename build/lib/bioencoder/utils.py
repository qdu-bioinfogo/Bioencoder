import os,time,sys
def getdatadir():
    basepath = os.path.abspath(__file__)
    folder = os.path.dirname(basepath)
    return folder


class FastaReader():
    file = None
    def __init__(self,file) -> None:
        self.file = file
    

    def readAll(self):     
        '''
        @msg: 读取一个fasta文件
        @param fa {str}  fasta 文件路径
        @return: {generator} 返回一个生成器，能迭代得到fasta文件的每一个序列名和序列
        '''
        with open(self.file,'r') as FA:
            seqName,seq='',''
            while 1:
                line=FA.readline()
                line=line.strip('\n')
                if (line.startswith('>') or not line) and seqName:
                    yield((seqName,seq))
                if line.startswith('>'):
                    seqName = line[1:]
                    seq=''
                else:
                    seq+=line
                if not line:break

    def getSeq(self,querySeqName,start=1,end=0):
        '''
        @msg: 获取fasta文件的某一条序列
        @param fa {str}  fasta 文件路径
        @param querySeqName {str}  序列名
        @param start {int}  截取该序列时，起始位置，可省略，默认为1
        @param end {int}  fasta 截取该序列时，最后位置，可省略，默认为该序列全长
        @return: {str} 返回找到(截取到)的序列
        '''
        if start<0: start=start+1
        for seqName,seq in self.readFa(self.file):
            if querySeqName==seqName:
                if end!=0: returnSeq = seq[start-1:end];print(start-1)
                else: returnSeq = seq[start-1:]
                return returnSeq

    def getReverseComplement(sequence):
        '''
        @msg: 获取DNA反向互补序列
        @param sequence {str}  一段DNA序列
        @return: {str} 返回反向互补序列
        '''
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()[::-1]

    def getGC(sequence):
        '''
        @msg: 获取某一条序列的GC含量
        @param sequence {str}  一段DNA序列
        @return: {float} 返回GC含量
        '''
        sequence=sequence.upper()
        return (sequence.count("G")+sequence.count("C"))/len(sequence)


    def getGapPos(sequence):
        '''
        @msg: 获取某条序列中gap的位置
        @param sequence {str}  一段DNA序列
        @return: {list}  返回一个列表，列表中每个元素为每个gap的起始和结束位置
        '''
        Ns = {'N', 'n'}
        result = []
        i = 0
        for base in sequence:
            i += 1
            if not base in Ns: continue
            if len(result) == 0 : result.append([i,i])
            elif i - result[-1][1] == 1: result[-1][1] = i
            else: result.append([i,i])
        return result

    def readSeqByWindow(sequence,winSize,stepSize):
        '''
        @msg: 滑窗读取某一条序列
        @param sequence {str}  一段DNA序列
        @param winSize {int}  窗口大小
        @param stepSize {int}  步长
        @return: {generator}  返回一个生成器，可迭代得到该序列的每一个窗口序列
        '''
        if stepSize<=0: return False
        now = 0
        seqLen = len(sequence)
        while(now+winSize-stepSize<seqLen):
            yield sequence[now:now+winSize]
            now+=stepSize
    

class ShowProcess():
    """
    显示处理进度的类
    调用该类相关函数即可实现处理进度的显示
    """
    i = 0 # 当前的处理进度
    max_steps = 0 # 总共需要处理的次数
    max_arrow = 50 #进度条的长度
    infoDone = 'done'

    # 初始化函数，需要知道总共的处理次数
    def __init__(self, max_steps, infoDone = 'Done'):
        self.max_steps = max_steps
        self.i = 0
        self.infoDone = infoDone

    # 显示函数，根据当前的处理进度i显示进度
    # 效果为[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]100.00%
    def show_process(self, i=None):
        if i is not None:
            self.i = i
        else:
            self.i += 1
        num_arrow = int(self.i * self.max_arrow / self.max_steps) #计算显示多少个'>'
        num_line = self.max_arrow - num_arrow #计算显示多少个'-'
        percent = self.i * 100.0 / self.max_steps #计算完成进度，格式为xx.xx%
        process_bar = '[' + '>' * num_arrow + '-' * num_line + ']'\
                      + '%.2f' % percent + '%' + '\r' #带输出的字符串，'\r'表示不换行回到最左边
        sys.stdout.write(process_bar) #这两句打印字符到终端
        sys.stdout.flush()
        if self.i >= self.max_steps:
            self.close()

    def close(self):
        print('')
        print(self.infoDone)
        self.i = 0

def Rvalue(aa1, aa2, AADict, Matrix):
	return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

class DimCounter():
    def EAAC(seqLength,winsize):
        return (seqLength - winsize + 1)*20


