#include<stdio.h>
#include<iostream>
#include<math.h>
#include<Windows.h>
#define pi 3.1415926
#define MaxSize 256
using namespace std;
struct HuffmanNode
{
	int weight;
	unsigned char data;
	HuffmanNode* llink;
	HuffmanNode* rlink;
};
struct HuffmanCode
{
	char bit[256]; //对于灰度图像，最多256个统计值，最多有256的编码   
	int len;	// 编码的实际长度
	unsigned char c;  // 原始字符
};
struct Stack
{
	char Data[MaxSize];
	int index;

	Stack() { this->index = -1; }

	void clear() { this->index = -1; }
	bool isFull() { return this->index >= MaxSize - 1; }
	bool isEmpty() { return this->index == -1; }
	void push(int d) {
		this->Data[++this->index] = d;
	}
	int pop() {
		return this->Data[this->index--];
	}
	int size() { return this->index + 1; }
};
struct MinHeap {

	HuffmanNode* Data[MaxSize];
	int Size;

	bool IsFull()
	{
		return this->Size == MaxSize - 1;
	}

	bool IsEmpty()
	{
		return this->Size == 0;
	}

	bool Push(HuffmanNode* X)
	{
		if (IsFull()) return false;

		int i;
		++this->Size;
		i = this->Size - 1;  //i指向插入后堆中的最后一个元素的位置
		for (; this->Data[(i - 1) / 2]->weight > X->weight && i > 0; i = (i - 1) / 2) //与父结点比较
		{
			this->Data[i] = this->Data[(i - 1) / 2];  //若大于父结点则与父结点交换
		}
		this->Data[i] = X;  //找到了合适的位置，将X插入

		return true;
	}

	void PercDown(int p)
	{
		HuffmanNode* X = this->Data[p]; //取出根结点的值

		int parent, child;
		for (parent = p; parent * 2 + 1 < this->Size; parent = child) { //当左子结点存在时继续
			child = parent * 2 + 1;
			if (child != this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // 右子结点存在，且小于左子结点
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // 找到了合适的位置
			else
				this->Data[parent] = this->Data[child]; // 下沉X
		}
		this->Data[parent] = X;
	}

	void BuildHeap() {
		for (int i = this->Size / 2 - 1; i >= 0; --i) PercDown(i);
	}

	HuffmanNode* Pop()
	{
		HuffmanNode* maxItem = this->Data[0]; //取出根结点存放的最大值
		HuffmanNode* X = this->Data[this->Size - 1];
		this->Size--;

		int parent, child;
		for (parent = 0; parent * 2 + 1 < this->Size; parent = child) { //当左子结点存在时继续
			child = parent * 2 + 1;
			if (child < this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // 右子结点存在，且大于左子结点
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // 找到了合适的位置
			else
				this->Data[parent] = this->Data[child];  // 下沉X   
		}
		this->Data[parent] = X;

		return maxItem;
	}
};
void deleteHuffmanTree(HuffmanNode* T);
HuffmanNode* createHuffmanTree(unsigned char buff[256], int weightArray[256], int n);
//创建huffman树，weightArray输入权重数组，n是待编码的序列长度
Stack codeStack;
int huffmanCoding(HuffmanNode* T, HuffmanCode* codeTable, int& index);
void calculateGrayLevelWeight(unsigned char* imgBuf, int width, int height, unsigned char grayLevel[],
	int grayLevelWeight[], int* grayLevelCount);
void generateFile(unsigned char* imgBuf, int width, int height, const char* fileout);//哈夫曼编码并写成文件
int decode8BitBuf(char codebuf, int ncode, HuffmanCode* codeTable, int matchPos[],
	unsigned char decodedBuf[], int decodedNum, int imgSize);
unsigned char* decoding(const char* filename, int* width, int* height);
bool saveBmp(const char* bmpName, unsigned char* imgBuf, int width, int height, int byteCount);
struct point//坐标点
{
	int x;
	int y;
};
void CircularFilter(int input[], double output[], int m, int n, 
	int filter_r, int N, double filter_thresh);//对灰度图像进行圆周滤波
int find_seed(double filter_out[], int m, int n, point p[],
	int  cluster_r, int  cluster_r_y);//寻找种子点
void regionGrowing_raw(unsigned char* buf, int m,
	int n, int seed_y, int seed_x, int thresh, unsigned char* bufOut);//区域生长法
unsigned char* readBmp(const char* bmpName, int* width, int* height, int* byteCount);//读入bmp文件




int main()
{
	//参数的设定，用于007
	int m = 985, n = 769;// 图片的长宽为m=985,n=769
	int filter_r = 20, N = 60;//圆周滤波的半径以及点数
	double filter_thresh = 0.4;//圆周滤波的阈值
	int cluster_r_x = 4, cluster_r_y = 4;//选取种子点预定区域的x、y宽度
	int thresh = 70;//种子生长的阈值


	//进行文件的读取
	FILE* fp;
	errno_t err = fopen_s(&fp, "data\\gray.raw", "r");//文件事先用处理为灰度图像并用PS保存成raw格式
	if (err != 0)
	{//检查文件是否打开
		printf("siu_fail");
		return 0;
	}
	unsigned char* bufGray = new unsigned char[m * n];//用于保存灰度图像
	fread(bufGray, sizeof(unsigned char), m * n, fp);//将图像读入bufGray中
	fclose(fp);
	int* buf = new int[m * n];//将无符号字符型转化为int类型的矩阵
	double* filter_result = new double [m*n];//圆周滤波后的矩阵
	for (int i = 0; i < m * n; i++)
	{//进行矩阵类型的转化已经filter_result的初始化
		buf[i] = (int)bufGray[i];
		filter_result[i] = 0.0;
	}
	CircularFilter(buf, filter_result,m,n,filter_r,N,filter_thresh);//进行圆周滤波
	printf("圆周滤波完成\n");
	point* position = new point[10000];//申请种子坐标矩阵
	int seed_num=find_seed(filter_result, m, n, position, cluster_r_x,cluster_r_y);//寻找种子生长点
	unsigned char* bufOut = new unsigned char[m * n];//输出的数组
	for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				bufOut[i * m + j] = 255;
			}
		}
	for (int i = 0; i < seed_num; i ++)
	{//对寻找到的种子点进行种子生长
		regionGrowing_raw(bufGray, m, n, position[i].x+3, 
			position[i].y-18, thresh, bufOut);
	}
	printf("种子生长完成\n");
	//err = fopen_s(&fp, "output.raw", "wb");
	//fwrite(bufOut, sizeof(unsigned char), m * n, fp);
	//fclose(fp);
	char fn[255] = "data\\outcome.huff";
	generateFile(bufOut, m, n, fn);//将最终的结果写入huff文件并保存
	printf("已经写入huff文件\n");

	/*char fy[255] = "data\\find_007.bmp";
	saveBmp(fy, bufOut, m, n, 1);
	printf("已经保存为bmp文件\n");*/
	/*下面进行huff文件的解码写入bmp文件*/
	printf("是否对huff文件进行解码？(Y/N):\n");
	char siusiu;
	cin >> siusiu;
	if (siusiu == 'Y')
	{//如果需要进行解码写成bmp
		char fx[255] = "data\\recover.bmp";
		int width, height, byteCount = 1;
		//解码
		unsigned char* buf2 = decoding(fn, &width, &height);
		printf("解码完成\n");
		saveBmp(fx, buf2, width, height, byteCount);
		printf("已经写为bmp\n");
	}

	//下面将滤波结果写入txt文件中
	
	//err = fopen_s(&fp, "012.txt", "w");//下面将滤波结果写入txt文件中
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < m; j++)
	//	{
	//		fprintf(fp, "%d%c", int(bufOut[m * i + j]), ' ');
	//	}
	//	fprintf(fp, "%c",'\n');
	//}
	//fclose(fp);
	
	delete[]bufGray;
	delete[]buf;
	delete[]filter_result;
	delete[]position;
	delete[]bufOut;
	
}




void CircularFilter(int input[], double output[], int m,
	int n,int filter_r, int N, double filter_thresh)
{	//参数依次为输入的图像、输出的图像、图像长、
	//图像宽、圆周滤波半径、选取的点数、滤波阈值
	long double A, B;
	double max = -1;//求取最大值，方便归一化
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			A = 0,B = 0;
			if (i - filter_r >= 0 && i + filter_r < n &&
				j - filter_r >= 0 && j + filter_r < m)
			{//判断是否在图像内
				for (int k = 0; k < N; k++)
				{//进行圆周滤波
					A = A + input[m*(int(i + filter_r * sin(2 * pi * k / N)))+
						int(j + filter_r * cos(2 * pi * k / N))] * cos(8 * pi * k / N);
					B = B + input[m*(int(i + filter_r * sin(2 * pi * k / N)))+
						int(j + filter_r * cos(2 * pi * k / N))] * sin(8 * pi * k / N);
				}
				output[i * m + j] = pow(A,2) + pow(B, 2);
				if (output[i * m + j] > max)//对图像求取最大值
					max = output[i * m+ j];
			}
		}
	}
	for (int i = 0; i < m * n; i++)
	{//对滤波后的图像进行二值化
		output[i] = output[i] / max;//对图像进行归一化
		if (output[i] < filter_thresh)//小于预订的阈值时为黑
			output[i] = 0;
		else
			output[i] = 255;
	}
}
int find_seed(double filter_out[], int m,int n,point p[], 
	int  cluster_r_x,int  cluster_r_y)
{	//寻找种子点,参数依次为输入的滤波后矩阵、
	//矩阵长、宽、输出的x坐标矩阵、y坐标矩阵、聚类的半径
	int top = -1;
	for (int i = cluster_r_x; i < n-cluster_r_x; i++)
	{
		for (int j = cluster_r_y; j < m-cluster_r_y; j++)
		{
			if (filter_out[(i - cluster_r_x)* m + j] != 0
				&& filter_out[(i + cluster_r_x) * m + j] != 0
				&& filter_out[i * m  + j-cluster_r_y] != 0 
				&& filter_out[i * m + j + cluster_r_y] != 0)
			{
				top++;
				p[top].x = i;
				p[top].y = j;
			}
		}
	}
	return top + 1;//返回种子点数组的大小
}
void regionGrowing_raw(unsigned char* buf, int m, int n, 
	int seed_y, int seed_x, int thresh, unsigned char* bufOut)
{//改造区域生长法，对无符号字符型数组进行生长,参数依次为输入的图像，
 //图像长、图像宽、种子点的xy坐标、生长阈值、输出矩阵
	int i, j, k;//循环变量
	//二维数组direction代表中心点8邻域坐标与该点在x和y方向上的偏移,
	//其中第一维为x方向的偏移,第二维为y方向的偏移
	int direction[8][2] = { { 0,1 },{ 1,1 },{ 1,0 },{ 1,-1 },
		{ 0,-1 },{ -1,-1 },{ -1,0 },{ -1,1 } };

	//栈申请，并定义栈顶指针
	point* stack = new point[m * n];
	int top = -1;//栈顶指针

	//当前正处理的点和弹出的点
	int currentPoint_x, currentPoint_y, popPoint_x, popPoint_y;

	int temp1, temp2;//临时变量,存放种子点以及待生长点的值

	//记录种子点的值
	temp1 = *(buf + seed_y * m + seed_x);

	//将给定种子点置标记0,入栈
	bufOut[seed_y * m + seed_x] = 0;
	top++;
	stack[top].x = seed_x;
	stack[top].y = seed_y;

	//堆栈
	while (top > -1) {
		//弹出栈顶元素,该元素已经生长过
		popPoint_x = stack[top].x;
		popPoint_y = stack[top].y;
		top--;

		//考察被弹出元素周围是否有没有生长的元素
		for (k = 0; k < 8; k++) {
			currentPoint_x = popPoint_x + direction[k][0];
			currentPoint_y = popPoint_y + direction[k][1];
			if (currentPoint_x >= m || currentPoint_x < 0 || currentPoint_y >= n || currentPoint_y < 0)
			{
				continue;
			}
			temp2 = buf[currentPoint_y * m + currentPoint_x];
			if (abs(temp1 - temp2) < thresh && bufOut[currentPoint_y * m + currentPoint_x] == 255)
			{
				top++;
				stack[top].x = currentPoint_x;
				stack[top].y = currentPoint_y;
				bufOut[currentPoint_y * m + currentPoint_x] = 0;
			}

		}//k

	}//while(top)

	 //清除缓冲区
	delete[]stack;
}
void deleteHuffmanTree(HuffmanNode* T)
{
	if (T == 0)//空树，返回
		return;
	struct HuffmanNode* s[200], * t;
	int top;
	top = 0;
	s[top] = T;//根结点进栈
	while (top > -1) {
		t = s[top];
		top--;
		if (t->rlink != 0) {//右子树进栈
			top++;
			s[top] = t->rlink;
		}
		if (t->llink != 0) {//左子树进栈
			top++;
			s[top] = t->llink;
		}
		delete t;//释放结点
	}
}
HuffmanNode* createHuffmanTree(unsigned char buff[256], int weightArray[256], int n)
{
	MinHeap heap;
	heap.Size = n;
	for (int i = 0; i < n; ++i)
	{
		heap.Data[i] = new HuffmanNode;
		heap.Data[i]->data = buff[i];
		heap.Data[i]->llink = nullptr;
		heap.Data[i]->rlink = nullptr;
		heap.Data[i]->weight = weightArray[i];
	}
	heap.BuildHeap();

	for (int i = 0; i < n - 1; ++i)
	{
		HuffmanNode* node = new HuffmanNode;
		node->llink = heap.Pop();
		node->rlink = heap.Pop();
		node->weight = node->llink->weight + node->rlink->weight;
		heap.Push(node);
	}
	return heap.Pop();
}
bool saveBmp(const char* bmpName, unsigned char* imgBuf, int width, int height, int byteCount)
{
	if (!imgBuf)
		return 0;

	//灰度图像颜色表空间1024，彩色图像没有颜色表
	int colorTable = 0;
	if (byteCount == 1) colorTable = 1024;

	//一行象素字节数为4的倍数
	int lineByte = (width * byteCount + 3) / 4 * 4;

	FILE* fp;
	fopen_s(&fp, bmpName, "wb");
	if (fp == 0) return 0;

	//填写文件头
	BITMAPFILEHEADER fileHead;
	fileHead.bfType = 0x4D42;
	fileHead.bfSize =
		sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colorTable + lineByte * height;
	fileHead.bfReserved1 = 0;
	fileHead.bfReserved2 = 0;
	fileHead.bfOffBits = 54 + colorTable;
	fwrite(&fileHead, sizeof(BITMAPFILEHEADER), 1, fp);

	// 填写信息头
	BITMAPINFOHEADER head;
	head.biBitCount = byteCount * 8;
	head.biClrImportant = 0;
	head.biClrUsed = 0;
	head.biCompression = 0;
	head.biHeight = height;
	head.biPlanes = 1;
	head.biSize = 40;
	head.biSizeImage = lineByte * height;
	head.biWidth = width;
	head.biXPelsPerMeter = 0;
	head.biYPelsPerMeter = 0;
	fwrite(&head, sizeof(BITMAPINFOHEADER), 1, fp);

	//颜色表拷贝  
	if (colorTable == 1024)
	{
		unsigned char table[1024];
		for (int i = 0; i < 256; i++)
		{
			*(table + i * 4 + 0) = i;
			*(table + i * 4 + 1) = i;
			*(table + i * 4 + 2) = i;
			*(table + i * 4 + 3) = 0;
		}
		fwrite(table, 1024, 1, fp);
	}

	//准备数据并写文件
	unsigned char* buf = new unsigned char[height * lineByte];
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * byteCount; j++)
			*(buf + i * lineByte + j) = *(imgBuf + i * width * byteCount + j);
	}
	fwrite(buf, height * lineByte, 1, fp);

	delete[]buf;

	fclose(fp);

	return 1;
}
int decode8BitBuf(char codebuf, int ncode, HuffmanCode* codeTable, int matchPos[],
	unsigned char decodedBuf[], int decodedNum, int imgSize)
{
	int i, j;
	int t, k, l;//
	for (i = 0; i < 8; i++)   //依次获得每个比特
	{
		t = (codebuf & 1 << i) != 0;//获得当前比特

		for (j = 0; j < ncode; j++) {//码表中第j个字符编码项
			if (matchPos[j] == -1)//已经失败，不再理会
				continue;

			if (codeTable[j].bit[matchPos[j]] == t) {//码表中待匹配的位置与当前解码相同				
				matchPos[j]++;
			}
			else {//没匹配上,置-1，继续下一个
				matchPos[j] = -1;
				continue;
			}

			if (matchPos[j] == codeTable[j].len) {//配上一个字符,一个数的码长不超过过码表数组的长度
				decodedBuf[decodedNum] = codeTable[j].c;
				decodedNum++;

				//matchPos要重新初始化，所有码重新匹配
				for (l = 0; l < ncode; l++)
					matchPos[l] = 0;

				break;//当匹配上一个字符时，码表内其他字符就不再匹配，跳出循环
			}
		}//j

		if (decodedNum == imgSize)
			break;
	}//i
	return decodedNum;//当前已解码个数
}
unsigned char* decoding(const char* filename, int* width, int* height)
{
	// 作业
	FILE* fp;
	fopen_s(&fp, filename, "rb");
	//fp = fopen(filename, "rb");
	int ncode;//编码表长
	fread(&ncode, sizeof(int), 1, fp);
	HuffmanCode huffmanCode[256];
	for (int i = 0; i < ncode; i++)//读出Huffman表
	{
		fread(&huffmanCode[i].c, sizeof(huffmanCode[i].c), 1, fp);
		fread(&huffmanCode[i].len, sizeof(huffmanCode[i].len), 1, fp);
		fread(&huffmanCode[i].bit, sizeof(char), huffmanCode[i].len, fp);
	}
	fread(width, sizeof(int), 1, fp);
	fread(height, sizeof(int), 1, fp);
	int w = *width;
	int h = *height;
	unsigned char* img_buf = new unsigned char[h * w];
	int matchPos[256] = { 0 };
	int decodedNum = 0;
	while (decodedNum < h * w)//依次解码每一个Huffman编码并存入buf中
		//每次读8位传给解码函数
	{
		char readbit;
		fread(&readbit, sizeof(char), 1, fp);
		//			printf("readbit: %d ",readbit);
		decodedNum = decode8BitBuf(readbit, ncode, huffmanCode, matchPos, img_buf, decodedNum, (*width) * (*height));
		//			printf("decodenum:%d\n",decodedNum);
	}

	//	return nullptr;
	return img_buf;
}
void generateFile(unsigned char* imgBuf, int width, int height, const char* fileout)
{
	//统计图像中频率非0的灰度级及其权重
	unsigned char grayLevel[256];
	int grayLevelWeight[256], grayLevelCount;
	calculateGrayLevelWeight(imgBuf, width, height, grayLevel, grayLevelWeight, &grayLevelCount);

	//创建哈夫曼树并生成编码数组
	HuffmanNode* HuffTree = createHuffmanTree(grayLevel, grayLevelWeight, grayLevelCount);//创建huffman树
	HuffmanCode huffmanCode[256];//存放哈夫曼编码
	codeStack.clear();
	int index = 0;
	int ncode = huffmanCoding(HuffTree, huffmanCode, index);//编码n个叶结点，huffmanTree树上一共有2n-1个点

	//创建文件
	FILE* fp;
	fopen_s(&fp, fileout, "wb");

	fwrite(&ncode, sizeof(ncode), 1, fp);
	for (int i = 0; i < ncode; ++i)
	{
		fwrite(&huffmanCode[i].c, sizeof(huffmanCode[i].c), 1, fp);
		fwrite(&huffmanCode[i].len, sizeof(huffmanCode[i].len), 1, fp);
		fwrite(&huffmanCode[i].bit, sizeof(char), huffmanCode[i].len, fp);
	}

	//图像宽、高写进文件，解码时才能申请合适的数据空间
	fwrite(&width, sizeof(int), 1, fp);
	fwrite(&height, sizeof(int), 1, fp);

	//对每一个字符编码写入的过程
	int i, j, k, t1, t2;
	char codebuf;//编码缓冲区
	t1 = 0;//缓冲区内编码位置
	for (i = height - 1; i >= 0; i--) {
		for (j = 0; j < width; j++) {
			//找到对应灰度级和权重
			for (k = 0; k < ncode; k++) {
				if (huffmanCode[k].c == imgBuf[i * width + j])
					break;
			}
			for (t2 = 0; t2 < huffmanCode[k].len; t2++) {//编码序列依次写入缓冲区
				if (huffmanCode[k].bit[t2] == 0)
					codebuf &= ~(1 << t1);//编码位置置0
				else
					codebuf |= (1 << t1);//编码位置置1
				t1++;
				if (t1 == 8) {//编码缓冲区满，写入文件
					fwrite(&codebuf, sizeof(char), 1, fp);
					t1 = 0;
				}
			}
			if (i == height - 1 && j == width - 1) {//到了最后一个像素
				if (t1 != 0)//最后一次缓冲区还没填满，写入文件，结束循环
					fwrite(&codebuf, sizeof(char), 1, fp);
			}
		}//j
	}//i
	fclose(fp);

	deleteHuffmanTree(HuffTree);
}
int huffmanCoding(HuffmanNode* T, HuffmanCode* codeTable, int& index)
{
	if (T->llink && T->rlink)  //Huffman树不存在度为1的结点
	{
		codeStack.push(0);
		huffmanCoding(T->llink, codeTable, index);
		codeStack.push(1);
		huffmanCoding(T->rlink, codeTable, index);
	}
	else
	{
		codeTable[index].len = codeStack.size();
		codeTable[index].c = T->data;
		for (int i = 0; i <= codeStack.index; ++i)
			codeTable[index].bit[i] = codeStack.Data[i];
		index++;
	}

	codeStack.pop();

	return index;	//返回code的数目
}
void calculateGrayLevelWeight(unsigned char* imgBuf, int width, int height, unsigned char grayLevel[],
	int grayLevelWeight[], int* grayLevelCount)
{
	int i, j, hist[256] = { 0 };
	for (i = 0; i < height; i++) {//统计直方图
		for (j = 0; j < width; j++) {
			hist[*(imgBuf + i * width + j)]++;
		}
	}
	for (i = 0; i < 256; i++)//统计频率/权重
		grayLevelWeight[i] = 0;
	int count = 0;
	for (i = 0; i < 256; i++) {//统计频率不为0的灰度级，用来生成哈夫曼树
		if (hist[i] != 0) {
			grayLevel[count] = i;
			grayLevelWeight[count] = hist[i];
			count++;
		}
	}
	*grayLevelCount = count;//实际灰度级数目
}
unsigned char* readBmp(const char* bmpName, int* width, int* height, int* byteCount)
{
	FILE* fp;
	fopen_s(&fp, bmpName, "rb");
	if (fp == 0) return 0;
	fseek(fp, sizeof(BITMAPFILEHEADER), 0);

	int w, h, b;
	BITMAPINFOHEADER head;
	fread(&head, sizeof(BITMAPINFOHEADER), 1, fp);
	w = head.biWidth;
	h = head.biHeight;
	b = head.biBitCount / 8;
	int lineByte = (w * b + 3) / 4 * 4;

	if (b == 1)
		fseek(fp, 1024, 1);
	unsigned char* imgBuf = new unsigned char[w * h * b];
	for (int i = 0; i < h; i++)
	{
		fread(imgBuf + i * w * b, w * b, 1, fp);
		fseek(fp, lineByte - w * b, 1);
	}
	fclose(fp);

	*width = w, * height = h, * byteCount = b;

	return imgBuf;
}