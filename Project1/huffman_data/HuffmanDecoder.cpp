#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#define nullptr NULL
#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),  (mode)))==NULL
#endif
#include<fstream>
#include<stddef.h>
struct HuffmanNode//Huffman���ڵ�
{
	int weight;
	unsigned char data;
	HuffmanNode *llink;
	HuffmanNode *rlink;
};

struct HuffmanCode
{
	char bit[256]; //���ڻҶ�ͼ�����256��ͳ��ֵ�������256�ı���
	int len;	// �����ʵ�ʳ���
	unsigned char c;  // ԭʼ�ַ�
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
		i = this->Size - 1;  //iָ��������е����һ��Ԫ�ص�λ��
		for (; this->Data[(i - 1) / 2]->weight > X->weight && i > 0; i = (i - 1) / 2) //�븸���Ƚ�
		{
			this->Data[i] = this->Data[(i - 1) / 2];  //��С�ڸ�������븸��㽻�� //С�ڣ������ϸ�����pzb
		}
		this->Data[i] = X;  //�ҵ��˺��ʵ�λ�ã���X����

		return true;
	}

	void PercDown(int p)//����ĳһ�����ڵ�
	{
		HuffmanNode* X = this->Data[p]; //ȡ��������ֵ

		int parent, child;
		for (parent = p; parent * 2 + 1 < this->Size; parent = child) { //�����ӽ�����ʱ����
			child = parent * 2 + 1;
			if (child != this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // ���ӽ����ڣ���С�����ӽ��
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // �ҵ��˺��ʵ�λ��
			else
				this->Data[parent] = this->Data[child]; // �³�X
		}
		this->Data[parent] = X;
	}

	void BuildHeap() {
		for (int i = this->Size / 2 - 1; i >= 0; --i) PercDown(i);
	}

	HuffmanNode* Pop()
	{
		HuffmanNode* minItem = this->Data[0]; //ȡ��������ŵ���Сֵ
		HuffmanNode* X = this->Data[this->Size - 1];
		this->Size--;

		int parent, child;
		for (parent = 0; parent * 2 + 1 < this->Size; parent = child) { //�����ӽ�����ʱ����
			child = parent * 2 + 1;
			if (child < this->Size - 1 && this->Data[child]->weight > this->Data[child + 1]->weight) { // ���ӽ����ڣ��Ҵ������ӽ��
				child++;
			}
			if (X->weight <= this->Data[child]->weight)
				break; // �ҵ��˺��ʵ�λ��
			else
				this->Data[parent] = this->Data[child];  // �³�X
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
			s[++top] = t;	// ��Ч�� top += 1; s[top] = t;
			t = t->llink;
		}

		if (top > -1) {
			t = s[top--]; // ��Ч�� t = s[top]; top -= 1;
			t = t->rlink;
		}
	}
}

void deleteHuffmanTree(HuffmanNode* T)
{
	if (T == 0)//����������
		return;
	struct HuffmanNode* s[200], * t;
	int top;
	top = 0;
	s[top] = T;//������ջ
	while (top > -1) {
		t = s[top];
		top--;
		if (t->rlink != 0) {//��������ջ
			top++;
			s[top] = t->rlink;
		}
		if (t->llink != 0) {//��������ջ
			top++;
			s[top] = t->llink;
		}
		delete t;//�ͷŽ��
	}
}

//����huffman����weightArray����Ȩ�����飬n�Ǵ���������г���
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
	if (T->llink && T->rlink)  //Huffman�������ڶ�Ϊ1�Ľ��
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
			codeTable[index].bit[i] = codeStack.Data[i];//�ӵ�λ��ʼ���λ����
														//����Huffman�������Ǵ��ڲ��������������벻һ��0�ǵ�����
		index++;
	}

	codeStack.pop();

	return index;	//����code����Ŀ
}

//��Ӧÿ��Ȩֵ�ı�������ʾ��T���������ĸ���㣬huffmanCode��Ӧ�ı���
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

//������ͼ����벿��
//ͳ��ͼ����ʵ�ʵĻҶȼ������ӦȨ�أ�grayLevel���ʵ�ʷ���Ҷȼ�������ֵ,
//grayLevelWeight��Ӧ��Ƶ�ʻ�Ȩ�أ�realGrayNum��0�Ҷȸ�����Ҳ��ǰ�������鳤��
void calculateGrayLevelWeight(unsigned char *imgBuf, int width, int height, unsigned char grayLevel[],
					  int grayLevelWeight[], int *grayLevelCount)
{
	int i, j, hist[256] = { 0 };
	for (i = 0; i<height; i++) {//ͳ��ֱ��ͼ
		for (j = 0; j<width; j++) {
			hist[*(imgBuf + i*width + j)]++;
		}
	}
	for (i = 0; i < 256; i++)//ͳ��Ƶ��/Ȩ��
		grayLevelWeight[i] = 0;
	int count = 0;
	for (i = 0; i < 256; i++) {//ͳ��Ƶ�ʲ�Ϊ0�ĻҶȼ����������ɹ�������
		if (hist[i] != 0) {
			grayLevel[count] = i;
			grayLevelWeight[count] = hist[i];
			count++;
		}
	}
	*grayLevelCount = count;//ʵ�ʻҶȼ���Ŀ
}

//���������벢д���ļ�
void generateFile(unsigned char *imgBuf, int width, int height, const char *fileout)
{
	//ͳ��ͼ����Ƶ�ʷ�0�ĻҶȼ�����Ȩ��
	unsigned char grayLevel[256];
	int grayLevelWeight[256], grayLevelCount;
	calculateGrayLevelWeight(imgBuf, width, height, grayLevel, grayLevelWeight, &grayLevelCount);

	//�����������������ɱ�������
	HuffmanNode *HuffTree = createHuffmanTree(grayLevel, grayLevelWeight, grayLevelCount);//����huffman��
	HuffmanCode huffmanCode[256];//��Ź���������
	codeStack.clear();
	int index = 0;
	int ncode = huffmanCoding(HuffTree, huffmanCode, index);//����n��Ҷ��㣬huffmanTree����һ����2n-1����
	showHuffmanCode(huffmanCode, ncode);//չʾ�����
	//�����ļ�
	FILE * fp;
//	fopen_s(&fp, fileout, "wb");
	fp = fopen(fileout, "wb");

	fwrite(&ncode, sizeof(ncode), 1, fp);
	for (int i = 0; i < ncode; ++i)
	{//д���ļ�ͷ
		fwrite(&huffmanCode[i].c, sizeof(huffmanCode[i].c), 1, fp);
		fwrite(&huffmanCode[i].len, sizeof(huffmanCode[i].len), 1, fp);
		fwrite(&huffmanCode[i].bit, sizeof(char), huffmanCode[i].len, fp);
	}

	//ͼ�����д���ļ�������ʱ����������ʵ����ݿռ�
	fwrite(&width, sizeof(int), 1, fp);
	fwrite(&height, sizeof(int), 1, fp);

	//��ÿһ���ַ�����д��Ĺ���
	int i, j, k, t1, t2;
	char codebuf;//���뻺����
	t1 = 0;//�������ڱ���λ��
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			//�ҵ���Ӧ�Ҷȼ���Ȩ��
			for (k = 0; k < ncode; k++) {
				if (huffmanCode[k].c == imgBuf[i*width + j])
					break;
			}
			for (t2 = 0; t2 < huffmanCode[k].len; t2++) {//������������д�뻺����
				if (huffmanCode[k].bit[t2] == 0)
					codebuf &= ~(1 << t1);//����λ����0
				else
					codebuf |= (1 << t1);//����λ����1
				t1++;
				if (t1 == 8) {//���뻺��������д���ļ�
					fwrite(&codebuf, sizeof(char), 1, fp);
					t1 = 0;
				}
			}
			if (i == height - 1 && j == width - 1) {//�������һ������
				if (t1 != 0)//���һ�λ�������û������д���ļ�������ѭ��
					fwrite(&codebuf, sizeof(char), 1, fp);
			}
		}//j
	}//i
	fclose(fp);

	deleteHuffmanTree(HuffTree);
}

/*
*  ��һ���ַ����������н���
*
*  ���룺
*  codebuf���ȴ�����Ļ�����
*  ncode�������������
*  huffmanCode�����������
*  macthPos����huffmanCode�ȳ������ƥ�䵽��λ�ã�-1ʱ����ûƥ���ϡ�Ӧ���ں������ʼ�����ҳ�ʼ����ֵ��ÿ��������ʼλ��0
*  decodedBuf��������������
*  decodedNum���ۻ�������ַ�����
*
*  �������ǰ�ѽ������
*/
int decode8BitBuf(char codebuf, int ncode, HuffmanCode *codeTable, int matchPos[],//matchPos 's size is 256
				  unsigned char decodedBuf[], int decodedNum, int imgSize)//decodedBuf ������ decodedNum ��������
{
	int i, j;
	int t, k, l;//
	for (i = 0; i < 8; i++)   //���λ��ÿ������
	{
		t = (codebuf & 1 << i) != 0;//��õ�ǰ����

		for (j = 0; j < ncode; j++) {//����е�j���ַ�������
			if (matchPos[j] == -1)//�Ѿ�ʧ�ܣ��������
				continue;

			if (codeTable[j].bit[matchPos[j]] == t) {//����д�ƥ���λ���뵱ǰ������ͬ
				matchPos[j]++;
			}
			else {//ûƥ����,��-1��������һ��
				matchPos[j] = -1;
				continue;
			}

			if (matchPos[j] == codeTable[j].len) {//����һ���ַ�,һ�������볤���������������ĳ���
				decodedBuf[decodedNum] = codeTable[j].c;
				decodedNum++;

				//matchPosҪ���³�ʼ��������������ƥ��
				for (l = 0; l < ncode; l++)
					matchPos[l] = 0;

				break;//��ƥ����һ���ַ�ʱ������������ַ��Ͳ���ƥ�䣬����ѭ��
			}
		}//j

		if (decodedNum == imgSize)
			break;
	}//i
	return decodedNum;//��ǰ�ѽ������
}


//����
unsigned char *decoding(const char * filename, int *width, int *height)
{
	// ��ҵ
	FILE* fp;
//	fopen_s(&fp, filename, "rb");
	fp = fopen(filename, "rb");
	int ncode;//�����
	fread(&ncode, sizeof(int), 1, fp);
	HuffmanCode huffmanCode[256];
	for (int i=0; i< ncode; i++)//����Huffman��
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
	while(decodedNum < h*w)//���ν���ÿһ��Huffman���벢����buf��
			//ÿ�ζ�8λ�������뺯��
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

	//�Ҷ�ͼ����ɫ��ռ�1024����ɫͼ��û����ɫ��
	int colorTable = 0;
	if (byteCount == 1) colorTable = 1024;

	//һ�������ֽ���Ϊ4�ı���
	int lineByte = (width * byteCount + 3) / 4 * 4;

	FILE* fp;
//	fopen_s(&fp, bmpName, "wb");
	fp = fopen(bmpName, "wb");
	if (fp == 0) return 0;

	//��д�ļ�ͷ
	BITMAPFILEHEADER fileHead;
	fileHead.bfType = 0x4D42;
	fileHead.bfSize =
		sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colorTable + lineByte * height;
	fileHead.bfReserved1 = 0;
	fileHead.bfReserved2 = 0;
	fileHead.bfOffBits = 54 + colorTable;
	fwrite(&fileHead, sizeof(BITMAPFILEHEADER), 1, fp);

	// ��д��Ϣͷ
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

	//��ɫ����
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

	//׼�����ݲ�д�ļ�
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

// ����һ��ͼ���ļ�����·��������ͼ�����ݡ�
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
    //����
    unsigned char *buf2 = decoding(fn, &width, &height);
    saveBmp(fx, buf2, width, height, byteCount);

	system("pause");

	return 0;
}
