#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#define nullptr NULL
#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),  (mode)))==NULL
#endif
#include<fstream>
#include<stddef.h>
struct HuffmanNode//Huffman树节点
{
	int weight;
	unsigned char data;
	HuffmanNode *llink;
	HuffmanNode *rlink;
};

struct HuffmanCode
{
	char bit[256]; //对于灰度图像，最多256个统计值，最多有256的编码
	int len;	// 编码的实际长度
	unsigned char c;  // 原始字符
};

#define MaxSize 256

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
			this->Data[i] = this->Data[(i - 1) / 2];  //若小于父结点则与父结点交换 //小于，大结点上浮——pzb
		}
		this->Data[i] = X;  //找到了合适的位置，将X插入

		return true;
	}

	void PercDown(int p)//调整某一子树节点
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
		HuffmanNode* minItem = this->Data[0]; //取出根结点存放的最小值
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

		return minItem;
	}
};

void preOrder_stack(HuffmanNode* T)
{
	if (T == 0)
		return;
	HuffmanNode* s[200];
	int top = -1;
	HuffmanNode* t = T;
	while (t != nullptr || top > -1)	//
	{
		while (t) {
			printf("%3d", t->weight);
			s[++top] = t;	// 等效于 top += 1; s[top] = t;
			t = t->llink;
		}

		if (top > -1) {
			t = s[top--]; // 等效于 t = s[top]; top -= 1;
			t = t->rlink;
		}
	}
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

//创建huffman树，weightArray输入权重数组，n是待编码的序列长度
HuffmanNode * createHuffmanTree(unsigned char buff[256], int weightArray[256], int n)
{
	MinHeap heap;
	heap.Size = n;
	for (int i = 0; i < n; ++i)
	{
		heap.Data[i] = new HuffmanNode;
		heap.Data[i]->data = buff[i];
		heap.Data[i]->llink = nullptr;
		heap.Data[i]->rlink = nullptr;
		heap.Data[i]->weight= weightArray[i];
	}
	heap.BuildHeap();

	for (int i = 0; i < n-1; ++i)
	{
		HuffmanNode *node = new HuffmanNode;
		node->llink = heap.Pop();
		node->rlink = heap.Pop();
		node->weight = node->llink->weight + node->rlink->weight;
		heap.Push(node);
	}
	return heap.Pop();
}

Stack codeStack;
int huffmanCoding(HuffmanNode *T, HuffmanCode *codeTable, int& index)
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
			codeTable[index].bit[i] = codeStack.Data[i];//从低位开始向高位读，
														//由于Huffman树可能是存在不定情况，因而编码不一定0是单独的
		index++;
	}

	codeStack.pop();

	return index;	//返回code的数目
}

//对应每个权值的编码结果显示，T哈夫曼树的根结点，huffmanCode对应的编码
void showHuffmanCode(HuffmanCode huffmanCode[], int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		printf("%d 's Huffman code is: ", huffmanCode[i].c);
		for (j = 0; j < huffmanCode[i].len; j++)
		{
			printf("%d", huffmanCode[i].bit[j]);//emmm
		}
		printf("\n");
	}
}

//以下是图像编码部分
//统计图像中实际的灰度级及其对应权重，grayLevel存放实际非零灰度级的索引值,
//grayLevelWeight对应的频率或权重，realGrayNum非0灰度个数，也是前边两数组长度
void calculateGrayLevelWeight(unsigned char *imgBuf, int width, int height, unsigned char grayLevel[],
					  int grayLevelWeight[], int *grayLevelCount)
{
	int i, j, hist[256] = { 0 };
	for (i = 0; i<height; i++) {//统计直方图
		for (j = 0; j<width; j++) {
			hist[*(imgBuf + i*width + j)]++;
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

//哈夫曼编码并写成文件
void generateFile(unsigned char *imgBuf, int width, int height, const char *fileout)
{
	//统计图像中频率非0的灰度级及其权重
	unsigned char grayLevel[256];
	int grayLevelWeight[256], grayLevelCount;
	calculateGrayLevelWeight(imgBuf, width, height, grayLevel, grayLevelWeight, &grayLevelCount);

	//创建哈夫曼树并生成编码数组
	HuffmanNode *HuffTree = createHuffmanTree(grayLevel, grayLevelWeight, grayLevelCount);//创建huffman树
	HuffmanCode huffmanCode[256];//存放哈夫曼编码
	codeStack.clear();
	int index = 0;
	int ncode = huffmanCoding(HuffTree, huffmanCode, index);//编码n个叶结点，huffmanTree树上一共有2n-1个点
	showHuffmanCode(huffmanCode, ncode);//展示编码表
	//创建文件
	FILE * fp;
//	fopen_s(&fp, fileout, "wb");
	fp = fopen(fileout, "wb");

	fwrite(&ncode, sizeof(ncode), 1, fp);
	for (int i = 0; i < ncode; ++i)
	{//写入文件头
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
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			//找到对应灰度级和权重
			for (k = 0; k < ncode; k++) {
				if (huffmanCode[k].c == imgBuf[i*width + j])
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

/*
*  对一个字符缓冲区进行解码
*
*  输入：
*  codebuf：等待解码的缓冲区
*  ncode：霍夫曼码表长度
*  huffmanCode：霍夫曼码表
*  macthPos：与huffmanCode等长，存放匹配到的位置，-1时代表没匹配上。应该在函数外初始化，且初始化的值是每个码表的起始位置0
*  decodedBuf：解码结果缓冲区
*  decodedNum：累积解码的字符个数
*
*  输出：当前已解码个数
*/
int decode8BitBuf(char codebuf, int ncode, HuffmanCode *codeTable, int matchPos[],//matchPos 's size is 256
				  unsigned char decodedBuf[], int decodedNum, int imgSize)//decodedBuf 解码结果 decodedNum 给个缓存
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


//解码
unsigned char *decoding(const char * filename, int *width, int *height)
{
	// 作业
	FILE* fp;
//	fopen_s(&fp, filename, "rb");
	fp = fopen(filename, "rb");
	int ncode;//编码表长
	fread(&ncode, sizeof(int), 1, fp);
	HuffmanCode huffmanCode[256];
	for (int i=0; i< ncode; i++)//读出Huffman表
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
	int matchPos[256]={0};
	int decodedNum = 0;
	while(decodedNum < h*w)//依次解码每一个Huffman编码并存入buf中
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
//	fopen_s(&fp, bmpName, "wb");
	fp = fopen(bmpName, "wb");
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

// 给定一个图像文件及其路径，读入图像数据。
unsigned char *readBmp(const char *bmpName, int *width, int *height, int *byteCount)
{
	FILE *fp;
//	fopen_s(&fp, bmpName, "rb");
	fp = fopen(bmpName, "rb");
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
	unsigned char *imgBuf = new unsigned char[w * h * b];
	for (int i = 0; i<h; i++)
	{
		fread(imgBuf + i*w*b, w*b, 1, fp);
		fseek(fp, lineByte - w*b, 1);
	}
	fclose(fp);

	*width = w, *height = h, *byteCount = b;

	return imgBuf;
}

int main()
{
	char fn[255]="outcome.huff";
	char fx[255]="recover.bmp";
    int width, height, byteCount=1;
    //解码
    unsigned char *buf2 = decoding(fn, &width, &height);
    saveBmp(fx, buf2, width, height, byteCount);

	system("pause");

	return 0;
}
