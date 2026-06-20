using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EmbeddingJPEG
{
    
    public class StegoImageJPEG
    {
        string filename;
        BinaryReader reader;
        int width = 0;
        int height = 0;
        int comp_count = 0;
        int[,] comp_params;
        int[,] mcu_params;
        Dictionary<string, byte>[] Decoding_Tables;
        Dictionary<byte, string>[] Encoding_Tables;
        int restartInterval = 0;
        int data_start_pos = 0;
        public short[,,] DCTCoeffs { get; private set; }

        public StegoImageJPEG(string filename)
        {
            this.filename = filename;
            reader = new BinaryReader(File.OpenRead(filename));
            Decoding_Tables = new Dictionary<string, byte>[4];
            for (int i = 0; i < 4; i++)
                Decoding_Tables[i] = new Dictionary<string, byte>();
        }

        public void Parse()
        {
            byte marker = 0x00;
            while (reader.BaseStream.Position < reader.BaseStream.Length)
            {
                marker = ReadMarker(reader);
                switch (marker)
                {
                    case 0xC0: // Сегмент начала кадра
                        ReadFrame();
                        break;
                    case 0xC2: // Отличительный маркер в файле формата Progressive JPEG
                        Console.WriteLine("Progressive JPEG format. Only baseline format supported.");
                        return;
                    case 0xC4: // Сегмент таблицы кодов Хаффмана
                        ReadHuffmanTable();
                        break;
                    case 0xDD: // Сегмент интервала перезапуска предсказания DC-коэффициентов
                        reader.ReadBytes(2);
                        restartInterval = (reader.ReadByte() << 8) | reader.ReadByte();
                        break;
                    case 0xDA: // Сегмент начала скана
                        ReadScan();
                        byte[] dataBytes = RemoveBitStaffing();
                        DCTCoeffs = ReadScanData(dataBytes);
                        break;
                    default: break;
                }
            }
        }
        /// <summary>
        /// Вспомогательный метод для чтения очередного маркера в бинарном файле
        /// </summary>
        /// <returns>Значение маркера без префикса</returns>
        private byte ReadMarker(BinaryReader reader)
        {
            byte prefix = reader.ReadByte();
            if (prefix == 0xFF)
                return reader.ReadByte();
            return 0x00;
        }
        /// <summary>
        /// Чтение сегмента начала кадра
        /// </summary>
        private void ReadFrame()
        {
            int length = (reader.ReadByte() << 8) | reader.ReadByte();
            int precision = reader.ReadByte();
            height = (reader.ReadByte() << 8) | reader.ReadByte();
            width = (reader.ReadByte() << 8) | reader.ReadByte();
            comp_count = reader.ReadByte();
            comp_params = new int[comp_count, 6];
            int compID;
            byte discretCoef;
            for (int i = 0; i < comp_count; i++)
            {
                compID = reader.ReadByte() - 1;
                discretCoef = reader.ReadByte();
                // Горизонтальный факторв
                comp_params[compID, 0] = (discretCoef >> 4) & 0x0F;
                // Вертикальный фактор
                comp_params[compID, 1] = discretCoef & 0x0F;
                // Число блоков данной компоненты
                comp_params[compID, 2] = comp_params[compID, 0] * comp_params[compID, 1];
                // Номер таблицы квантования
                comp_params[compID, 3] = reader.ReadByte();

            }
        }
        /// <summary>
        /// Считывание таблицы кодов Хаффмана
        /// </summary>
        private void ReadHuffmanTable()
        {
            int length = (reader.ReadByte() << 8) | reader.ReadByte() - 19;
            byte tableInfo = reader.ReadByte();
            // Номер пары таблиц (DC - AC)
            int tableNumber = tableInfo & 0x0F;
            // Тип таблицы: 0 = DC, 1 = AC
            int tableType = (tableInfo >> 4) & 0x0F;
            // Индекс в массиве таблиц
            int tableInd = (tableNumber * 2 + tableType);
            // Данные о количестве кодов каждой длины
            byte[] codeCounts = new byte[16];
            for (int i = 0; i < codeCounts.Count(); i++)
                codeCounts[i] = reader.ReadByte();
            // Значения для кодов
            byte[] values = new byte[length];
            for (int i = 0; i < values.Count(); i++)
                values[i] = reader.ReadByte();
            string codeString;
            string[] codes = new string[length];
            int code = 0;
            int k = 0;
            // Формирование кода по стандартному алгоритму
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < codeCounts[i]; j++)
                {
                    codeString = Convert.ToString(code, 2);
                    codes[k] = new string('0', i + 1 - codeString.Length) + codeString;
                    code++;
                    k++;
                }
                code = code << 1;
            }
            for (int i = 0; i < length; i++)
                Decoding_Tables[tableInd][codes[i]] = values[i];
        }
        /// <summary>
        /// Считывание начала сегмента скана
        /// </summary>
        private void ReadScan()
        {
            int length = (reader.ReadByte() << 8) | reader.ReadByte();
            int componentsCount = reader.ReadByte();
            int compID;
            byte TableNumbers;
            for (int i = 0; i < comp_count; i++)
            {
                compID = reader.ReadByte() - 1;
                TableNumbers = reader.ReadByte();
                // Индекс таблицы DC
                comp_params[compID, 4] = (TableNumbers >> 4) & 0x0F;
                // Индекс таблицы AC
                comp_params[compID, 5] = TableNumbers & 0x0F;
            }
            // Три байта, завершающие начало сегмента
            reader.ReadBytes(3);
        }
        /// <summary>
        /// Удаление стаффинга нулевых байт в битовом потоке данных
        /// </summary>
        /// <returns>Данные, собранные в массив байт</returns>
        private byte[] RemoveBitStaffing()
        {
            // Позиция начала данных в файле
            data_start_pos = (int)reader.BaseStream.Position;
            long length = reader.BaseStream.Length - reader.BaseStream.Position - 2;
            bool staffing = false;
            byte[] data = new byte[length];
            byte currentByte;
            int counter = 0;
            for (int i = 0; i < length; i++)
            {
                currentByte = reader.ReadByte();
                if (!staffing)
                {
                    if (currentByte == 0xFF)
                    {
                        data[counter] = currentByte;
                        staffing = true;
                    }
                    else
                    {
                        data[counter] = currentByte;
                    }
                    counter++;
                }
                else
                {
                    if (currentByte == 0x00)
                    {
                        staffing = false;
                    }
                    else if (currentByte >= 0xD0 && currentByte <= 0xD7) // Маркеры перезапуска
                    {
                        data[counter] = currentByte;
                        staffing = false;
                        counter++;
                    }
                    else
                    {
                        data[counter] = currentByte;
                        counter++;
                    }
                }
            }
            return data;
        }
        /// <summary>
        /// Считывание битового потока энтропийно закодированных данных скана
        /// </summary>
        /// <param name="dataBytes">Данные в виде массива байт</param>
        /// <returns>Массив коэффициентов дискретного косинусного преобразования</returns>
        private short[,,] ReadScanData(byte[] dataBytes)
        {
            int blocksPerMCU = Math.Max(comp_params[0, 2], Math.Max(comp_params[1, 2], comp_params[2, 2]));
            int MCUSize = comp_params[0, 2] + comp_params[1, 2] + comp_params[2, 2];
            mcu_params = new int[4, MCUSize];
            int blockIndex = -1;
            for (int i = 0; i < comp_count; i++)
            {
                for (int j = 0; j < comp_params[i, 2]; j++)
                {
                    blockIndex++;
                    // DC table index in table array
                    mcu_params[0, blockIndex] = 2 * comp_params[i, 4];
                    // AC table index in table array
                    mcu_params[1, blockIndex] = 2 * comp_params[i, 5] + 1;
                    // Block index in component
                    mcu_params[2, blockIndex] = j;
                    // for first block in component
                    if (j == 0)
                        // last index in this component
                        mcu_params[3, blockIndex] = blockIndex + comp_params[i, 2] - 1;
                }
            }
            int MCUCount = height * width / 64 / blocksPerMCU;
            short[,,] data = new short[MCUCount, MCUSize, 64];
            byte buffer;
            string strBuffer;
            string accumulator = "";
            bool AC = false;
            bool readValue = false;
            byte value;
            int valueBitsCount = 0;
            int num;
            int size = 0;
            int byteIndex = -1;
            int coefIndex = 0;
            blockIndex = 0;
            int mcuIndex = 0;
            bool restrartMarkers = false;
            int mcuTillRestart = -1;

            if (restartInterval > 0)
            {
                restrartMarkers = true;
                mcuTillRestart = restartInterval;
            }

            while (true)
            {
                byteIndex++;
                buffer = dataBytes[byteIndex];
                strBuffer = Convert.ToString(buffer, 2);
                strBuffer = new string('0', 8 - strBuffer.Length) + strBuffer;
                for (int j = 0; j < 8; j++)
                {
                    accumulator = accumulator + strBuffer.Substring(0, 1);
                    strBuffer = strBuffer.Substring(1, strBuffer.Length - 1);
                    if (!readValue)
                    {
                        if (AC)
                        {
                            if (Decoding_Tables[mcu_params[1, blockIndex]].ContainsKey(accumulator))
                            {
                                if (Decoding_Tables[mcu_params[1, blockIndex]][accumulator] == 0x00)
                                {
                                    coefIndex = 64;
                                    AC = false;
                                    accumulator = "";
                                }
                                else if (Decoding_Tables[mcu_params[1, blockIndex]][accumulator] == 0xF0)
                                {
                                    coefIndex += 16;
                                    accumulator = "";
                                }
                                else
                                {
                                    value = Decoding_Tables[mcu_params[1, blockIndex]][accumulator];
                                    coefIndex += value >> 4;
                                    size = value & 0x0F;
                                    valueBitsCount = size;
                                    accumulator = "";
                                    readValue = true;
                                }
                            }
                        }
                        else
                        {
                            if (Decoding_Tables[mcu_params[0, blockIndex]].ContainsKey(accumulator))
                            {
                                if (Decoding_Tables[mcu_params[0, blockIndex]][accumulator] == 0x00)
                                {
                                    if (mcuIndex > 0 && mcuTillRestart < restartInterval && mcu_params[2, blockIndex] == 0)
                                        data[mcuIndex, blockIndex, coefIndex] =
                                            data[mcuIndex - 1, mcu_params[3, blockIndex], coefIndex];
                                    else if (mcu_params[2, blockIndex] > 0)
                                        data[mcuIndex, blockIndex, coefIndex] =
                                            data[mcuIndex, blockIndex - 1, coefIndex];
                                    else
                                        data[mcuIndex, blockIndex, coefIndex] = 0;
                                    coefIndex++;
                                    accumulator = "";
                                    AC = true;
                                }
                                else
                                {
                                    size = Decoding_Tables[mcu_params[0, blockIndex]][accumulator];
                                    valueBitsCount = size;
                                    accumulator = "";
                                    readValue = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        if (AC)
                        {
                            valueBitsCount--;
                            if (valueBitsCount <= 0)
                            {
                                num = Convert.ToInt32(accumulator, 2);
                                data[mcuIndex, blockIndex, coefIndex] = (short)num;
                                if ((num >> (size - 1)) == 0)
                                    data[mcuIndex, blockIndex, coefIndex] -= (short)((1 << size) - 1);
                                accumulator = "";
                                readValue = false;
                                coefIndex++;
                            }
                        }
                        else
                        {
                            valueBitsCount--;
                            if (valueBitsCount <= 0)
                            {
                                num = Convert.ToInt32(accumulator, 2);
                                if (mcuIndex > 0 && mcuTillRestart < restartInterval && mcu_params[2, blockIndex] == 0)
                                    data[mcuIndex, blockIndex, coefIndex] =
                                        (short)(data[mcuIndex - 1, mcu_params[3, blockIndex], coefIndex] + num);
                                else if (mcu_params[2, blockIndex] > 0)
                                    data[mcuIndex, blockIndex, coefIndex] =
                                        (short)(data[mcuIndex, blockIndex - 1, coefIndex] + num);
                                else
                                    data[mcuIndex, blockIndex, coefIndex] = (short)num;
                                if ((num >> (size - 1)) == 0)
                                    data[mcuIndex, blockIndex, coefIndex] -= (short)((1 << size) - 1);
                                accumulator = "";
                                readValue = false;
                                coefIndex++;
                                AC = true;
                            }
                        }
                    }
                    if (coefIndex == 64)
                    {
                        blockIndex++;
                        if (blockIndex == MCUSize)
                        {
                            mcuIndex++;
                            blockIndex = 0;
                            if (restrartMarkers)
                            {
                                mcuTillRestart--;
                                if (mcuTillRestart == 0)
                                {
                                    strBuffer = "";
                                    byteIndex += 2;
                                    mcuTillRestart = restartInterval;
                                    coefIndex = 0;
                                    AC = false;
                                    break;
                                }
                            }
                        }
                        coefIndex = 0;
                        AC = false;

                    }
                    if (mcuIndex >= MCUCount)
                        return data;
                }
            }
        }
        /// <summary>
        /// Вспомогательный метод для инвертирования словарей - таблиц кодов Хаффмана
        /// </summary>
        private Dictionary<byte, string>[] ReverseHuffmanTables()
        {
            Dictionary<byte, string>[] Rev_Huffman_Tables = new Dictionary<byte, string>[Decoding_Tables.Length];
            for (int i = 0; i < Decoding_Tables.Length; i++)
            {
                Rev_Huffman_Tables[i] = new Dictionary<byte, string>();
                foreach (KeyValuePair<string, byte> pair in Decoding_Tables[i])
                    Rev_Huffman_Tables[i][pair.Value] = pair.Key;
            }
            return Rev_Huffman_Tables;
        }
        /// <summary>
        /// Кодирование и запись в файл массива коэффициентов
        /// дискретного косинусного преобразования
        /// </summary>
        /// <param name="data">Массив коэффициентов</param>
        /// <param name="newFilename">Имя нового файла для записи</param>
        public void WriteNewData(short[,,] data, string newFilename)
        {
            bool AC = false;
            byte size;
            int num;
            int run = 0;
            string numStr = "";
            string strData;
            int mcuTillRestart = -1;
            byte restartMarker = 0xD0;
            bool restartMarkers = false;
            string accumulator = "";
            if (restartInterval > 0)
            {
                restartMarkers = true;
                mcuTillRestart = restartInterval;
            }
            Encoding_Tables = ReverseHuffmanTables();
            File.Copy(filename, newFilename, overwrite: true);
            BinaryWriter writer = new BinaryWriter(File.OpenWrite(newFilename));
            writer.BaseStream.Seek(data_start_pos, SeekOrigin.Begin);
            for (int mcuIndex = 0; mcuIndex < data.GetLength(0); mcuIndex++)
            {
                for (int blockIndex = 0; blockIndex < data.GetLength(1); blockIndex++)
                {
                    for (int coeffIndex = 0; coeffIndex < data.GetLength(2); coeffIndex++)
                    {
                        if (AC)
                        {
                            num = data[mcuIndex, blockIndex, coeffIndex];
                            if (num == 0)
                            {
                                if (coeffIndex == 63)
                                {
                                    accumulator += Encoding_Tables[mcu_params[1, blockIndex]][0];
                                    run = 0;
                                }
                                else
                                    run++;
                            }
                            else
                            {
                                if (run > 15)
                                {
                                    while (run > 15)
                                    {
                                        accumulator += Encoding_Tables[mcu_params[1, blockIndex]][0xF0];
                                        run -= 16;
                                    }
                                }
                                size = (byte)(Math.Floor(Math.Log2(Math.Abs(num))) + 1);
                                if (num < 0)
                                    num += ((1 << size) - 1);
                                numStr = Convert.ToString(num, 2);
                                accumulator += Encoding_Tables[mcu_params[1, blockIndex]][(byte)((run << 4) | size)];
                                accumulator += new string('0', size - numStr.Length) + numStr;
                                run = 0;
                            }
                        }
                        else
                        {
                            if (mcuIndex > 0 && mcuTillRestart < restartInterval && mcu_params[2, blockIndex] == 0)
                                num = data[mcuIndex, blockIndex, coeffIndex] -
                                    data[mcuIndex - 1, mcu_params[3, blockIndex], coeffIndex];
                            else if (mcu_params[2, blockIndex] > 0)
                                num = data[mcuIndex, blockIndex, coeffIndex] -
                                    data[mcuIndex, blockIndex - 1, coeffIndex];
                            else
                                num = data[mcuIndex, blockIndex, coeffIndex];

                            if (num == 0)
                                accumulator += Encoding_Tables[mcu_params[0, blockIndex]][0];
                            else
                            {
                                size = (byte)(Math.Floor(Math.Log2(Math.Abs(num))) + 1);
                                if (num < 0)
                                    num += ((1 << size) - 1);
                                numStr = Convert.ToString(num, 2);
                                accumulator += Encoding_Tables[mcu_params[0, blockIndex]][size];
                                accumulator += new string('0', size - numStr.Length);
                                accumulator += numStr;
                            }
                            AC = true;
                        }
                        while (accumulator.Length > 7)
                        {
                            strData = accumulator.Substring(0, 8);
                            if (strData == "11111111")
                            {
                                writer.Write([0xFF]);
                                writer.Write([0x00]);
                            }
                            else
                            {
                                writer.Write(Convert.ToByte(strData, 2));
                            }
                            accumulator = accumulator.Substring(8);
                        }
                    }
                    AC = false;
                }
                if (restartMarkers && mcuIndex < data.GetLength(0) - 1)
                {
                    mcuTillRestart--;
                    if (mcuTillRestart == 0)
                    {
                        mcuTillRestart = restartInterval;
                        if (accumulator.Length > 0)
                        {
                            accumulator += new string('1', 8 - accumulator.Length);
                            writer.Write([Convert.ToByte(accumulator, 2)]);
                        }
                        writer.Write([0xFF]);
                        writer.Write([restartMarker]);
                        accumulator = "";
                        if (restartMarker < 0xD7)
                            restartMarker++;
                        else
                            restartMarker = 0xD0;
                    }
                }
                else if (mcuIndex == data.GetLength(0) - 1)
                {
                    if (accumulator.Length > 0)
                    {
                        accumulator += new string('1', 8 - accumulator.Length);
                        writer.Write(Convert.ToByte(accumulator, 2));
                    }
                    writer.Write([0xFF]);
                    writer.Write([0xD9]);
                }
            }
            writer.Close();
        }
    }
}
