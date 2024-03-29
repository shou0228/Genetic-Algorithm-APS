# Genetic-Algorithm for APS #自動化生產線排程 #平行機台 

✒️排程目標: 
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
隨著工單逐漸多樣化、產量少且交期短，多機同步作業與換線需求大大增加了排程的難度，使得原本的排程規劃已無法有效提升生產效率，甚至

產生閒置，因此我們打算以用 GA 演算法作為 APS (先進規劃排程)的演算邏輯，以 Python 程式語言進行實作來解決此排程難題。

✒️生產線製程資訊與訂單資訊:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
![image](https://user-images.githubusercontent.com/68886395/158193978-0402b276-8bfe-423b-9c65-15ba9304b01a.png)圖1 

![image](https://user-images.githubusercontent.com/68886395/158196237-71a49176-3093-4932-bbab-f3a46171610e.png)圖2

![image](https://user-images.githubusercontent.com/68886395/158204879-2d0fce92-8c37-4f38-82c0-6cb6fcc82f34.png)圖3

圖1高階顯卡種類三種，以 A、B、C 代表。各類顯卡工序相同（依序為自動鎖付風扇、人工鎖付背板、光學檢測步驟），但各工作站所需工作

時長不同且具有不同機種間換線時間。

圖2訂單資訊包含多筆工單，模擬生產現場之商品種類多、交期不同、數量不同。

圖3可得自動螺絲機一台人工鎖付背板之作業員數量為兩人，光學檢測 AOI 機兩台，性能及工作時長相同。

✒️ 排程問題描述:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
由於上述所提供生產線上製程相關資訊與訂單種類包含了換線時間、交期、工序與顯卡等諸多問題，我們為了讓此生產線達到生產時間最佳化

結合生產管理與作業研究相關知識，並結合個人的程式能力python語法以基因演算法撰寫此排程相關解決方式。

✒️ 程式導入前置:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>編碼原則:

![image](https://user-images.githubusercontent.com/68886395/158203159-c700fcb3-7e65-44f5-8325-04214134d3f6.png)

>生產線製程資訊轉換excel表(平行機台):

M0: 1台自動螺絲機、M1': 2位作業員 和 M2': 2台AOI機。各機台M數字部分為時間(分鐘)，Quantity 為各自工單的生產數量(個)

![image](https://user-images.githubusercontent.com/68886395/158203628-b995496e-0104-4e72-aff9-0763f71dc580.png)

該excel讀檔連結 [fsp.xlsx](https://github.com/shou0228/Genetic-Algorithm-APS/blob/e5435e57e00267109aecd5d84e571929c022ee4a/fsp.xlsx) 檔案內點選parallel

✒️ 基因演算法程式說明:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
a.初始設定:將讀檔後的資料及所需變數進行命名

![image](https://user-images.githubusercontent.com/68886395/158207172-2f2f9023-e779-4f31-98ef-2effae28b7fe.png)

“num_mc”為機器數，“num_job”為工件數，“pt”為處理時間的命名，“ct”為換線時間，主要是將資料內的位置以陣列方式去做讀取。

“typeSequence”主要為判斷工件 ABC 的種類。圖中 26-31 行為屆時程式可以自行輸入的觀察值、交配率、突變率、突變選擇率和所

需交配的次數等，也可選擇不輸入直接以本研究預設的默認值為準。

b.編碼:

編碼程序

![image](https://user-images.githubusercontent.com/68886395/158208361-cc726e95-33e4-49b7-8c91-f7ef77249038.png)

透過隨機的方式，產生初始序列的染色體。

![image](https://user-images.githubusercontent.com/68886395/158208411-14fb3725-e97d-4c46-9992-f8421d4bd714.png)

c.交配:

本研究藉由上一步驟產生的隨機序列染色體中，隨機抽取兩個兩個一組，根據初始設定的交配率來決定是否要進行交配，如果要交配

則藉由切割方式將兩組剛剛抽的染色體設為母代，並藉由切割染色體內基因線段的方式重新組成兩個子代。

![image](https://user-images.githubusercontent.com/68886395/158209023-3c7711cd-9d64-4a6e-893e-32e1519c503b.png)

![image](https://user-images.githubusercontent.com/68886395/158209048-d5e807aa-bac6-40e0-a831-835bcdd84809.png)

d.修復:

前一步驟所進行的交配，為隨機切割排列的方式。此會導致有些染色體內的基因(工件)出現次數過多或過少，如下圖所示，子代的5和3

會因此形成一個不可行的子代染色體排程解，此時本研究需將不可行的子代染色體進行修復動作，藉以得到可行解。

![image](https://user-images.githubusercontent.com/68886395/158209637-476e54e0-bb5b-4ffc-9567-a2fb229ce45a.png)

![image](https://user-images.githubusercontent.com/68886395/158209691-e8a94a38-a78d-4b44-83b2-2d9b06996166.png)

e.突變

依據先前突變選擇率(mutation selection rate)決定染色體中有多少百分比的基因要進行突變，假設每條染色體有 9 個基因，如下圖所示。

當 mutation selection rate 為 0.2 時，經由 9*0.2的計算並採四捨五入方式決定2個基因要進行突變。

![image](https://user-images.githubusercontent.com/68886395/158210506-10365c6d-5134-4067-980b-11909012b4a7.png)

![image](https://user-images.githubusercontent.com/68886395/158210538-edc1b6e1-20e2-4666-ae37-b72e14bd66be.png)

f.適應值計算

此問題目標為求出最短完工時間，故在適應值計算階段需寫出適合度函數的程式碼，用以計算適應值，即最短完工時間。

本階段程式碼分為 4 個部分，依序為「計算前設定」、「以機器做工件排序」、「排程邏輯」以及「適應值紀錄」。

(a) 計算前設定:

目的在於初始化適合度函數所需要的初始值以及串列儲存空間。其中下圖第 116 行的「total_chromsome」為儲存所有複

製的親代以及複製的子代染色體，用以被第 121 行的迴圈迭代，取出單一染色體用來被其他部分所計算。

![image](https://user-images.githubusercontent.com/68886395/158212624-661de217-6048-4ca6-8669-2c4e3509b5e4.png)

(b) 以機器做工件排序

此部分程式碼旨在處理以「機器分類」做工件排序，用以方便計算適應值

![image](https://user-images.githubusercontent.com/68886395/158213046-986466b6-b36d-412b-9995-c95a492d11d0.png)

0 表示為 Machine0 即自動鎖螺絲機，而值為一串列，”[2, 0, 1]”即為Machine0 的工件排序，其餘鍵值依此類推

(c) 排程邏輯

上圖第 137 行及第 138 行為本部分前置作業，第 137 行為抓出 Excel 報表中第 order 列，用以記錄各工單種類，以便在程式中判斷。

第 138 行為記錄每一個工單的作業後時間，時間以第一個工單為 0 算起。138-157 行為本演算法計算排程的核心部分「排程邏輯」，排程

邏輯分三個處理階段，分別是「若前一機器作業時間大於當前機器前一份工件的作業時間再加上換線時間，則延續前一機器作業時間。」、

「當前機器前一份工件加上當前工件，並確保不同工單之間會有換線時間」以及「若不滿足上述兩項，則前一機器作業時間直接加上當前機器作

業時間。」

(d) 適應值記錄

本部分程式碼主要記錄排程後的總完工時間（makespan），以利在後續的基因演算法「選擇」與「比較」階段取出使用。

g.選擇

本階段程式碼採用輪盤法(Roulette wheel)選擇較好的染色體，擁有較高適應值的染色體在輪盤中會有較大的比率被選擇

![image](https://user-images.githubusercontent.com/68886395/158215388-910b5a8d-c2b4-4bdb-b034-3a2b42722357.png)

h. 比較

先比較每條染色體的總加權延遲 (chrom_fit) ，選出此輪找到的最好解 (Tbest_now) ，接著在跟目前為止找到的最好解 (Tbest) 進行

比較，一旦這一輪的解比目前為止找到的解還要好，就替代 Tbest 並記錄該解所得到的排程結果，參考下圖為本階段程式碼終止條件。

![image](https://user-images.githubusercontent.com/68886395/158214132-411ec897-adae-4f45-9686-5205f3fe3bb2.png)

本研究使用 optimal_time_ratio 來當作本研究終止條件的方法，這方式主要使用當前找到的最佳解時間除以程式碼總經過執行時間，

若迭代小於 0.5 則終止迭代(iteration)，終止的近似最佳解即為結果。因為需要一定量的樣本，這裡的 n 設於大於 5000，

意味著染色體要交配超過 5000 次。

![image](https://user-images.githubusercontent.com/68886395/158214388-b999d3fe-eb7c-418d-a0dd-de768023ef4a.png)

i.結果

本結果列出(Print)五個內容，由上而下分別是工單排程（optimal sequence）、近似最佳解（optimal value）、總執行時間（the total time）、

近似最佳解時間（the optimal time）以及終止條件所使用的最佳解時間比率（the optimal time ratio）

參照下圖印出內容以及 XY關係圖（交配次數（generation）/總完工時間（makespan））的程式其中本 XY 關係圖表示交配50次數到 4000 餘次

即找到近似最佳解，而總交配次數為 8000 餘次，終止條件成立(optimal_time_ration < 0.5)，故交配跌代(iteration)停止。

![image](https://user-images.githubusercontent.com/68886395/158214947-d5b3d8b5-3ea9-405b-acd2-7c460ff438e2.png)

最終程式碼執行結果:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#9000/50000交配 ，可得最佳時間:958min

![image](https://user-images.githubusercontent.com/68886395/158214998-736d769d-9330-4615-a114-dfc951cea7a8.png)


















