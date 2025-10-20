package gobigwig

import (
	"bytes"
	"errors"
	"io"
	"net/http"
	"os"
	"strconv"
)

// bigWigFileType 表示文件类型
type bigWigFileType int

const (
	BWG_FILE bigWigFileType = iota
	BWG_HTTP
	BWG_HTTPS
	BWG_FTP
)
type size_t =int64


type URL struct {
	rs io.ReadSeeker // 实际用于 Read/Seek 的接口
	// 远程文件专用
	client *http.Client
	url    string
	buf    *bytes.Buffer
	Type         bigWigFileType
	FName        string
	IsCompressed bool
	FilePos      int64
}

// Open 打开本地文件或远程 URL
func Open(fname string) (*URL, error) {
	u := &URL{
		FName: fname,
	}
	switch {
	case len(fname) >= 7 && fname[:7] == "http://":
		u.Type = BWG_HTTP
		u.client = &http.Client{}
		u.url = fname
		u.buf = bytes.NewBuffer(nil)
		u.rs = u // 使用自定义 ReadSeeker
	case len(fname) >= 8 && fname[:8] == "https://":
		u.Type = BWG_HTTPS
		u.client = &http.Client{}
		u.url = fname
		u.buf = bytes.NewBuffer(nil)
		u.rs = u
	default:
		// 本地文件
		f, err := os.Open(fname)
		if err != nil {
			return nil, err
		}
		u.Type = BWG_FILE
		u.rs = f
	}

	return u, nil
}

// Close 关闭文件
func (u *URL) Close() error {
	if u.Type == BWG_FILE {
		if f, ok := u.rs.(*os.File); ok {
			return f.Close()
		}
	}
	// 远程文件没有长连接需要关闭
	return nil
}

// Read 实现 io.Reader
func (u *URL) Read(p []byte) (int, error) {
	if u.Type == BWG_FILE {
		return u.rs.Read(p)
	}
	// 远程文件，缓冲区读取
	if u.buf.Len() == 0 {
		err := u.fillBuffer()
		if err != nil {
			return 0, err
		}
	}
	return u.buf.Read(p)
}

// Seek 实现 io.Seeker
func (u *URL) Seek(offset int64, whence int) (int64, error) {
	if u.Type == BWG_FILE {
		return u.rs.Seek(offset, whence)
	}
	// 远程文件，通过 Range 请求实现
	var absPos int64
	switch whence {
	case io.SeekStart:
		absPos = offset
	case io.SeekCurrent:
		absPos = u.FilePos + offset
	case io.SeekEnd:
		// 远程文件不能直接 SeekEnd，需要提前知道长度
		return 0, errors.New("SeekEnd not supported for remote files")
	default:
		return 0, errors.New("invalid whence")
	}
	u.FilePos = absPos
	u.buf.Reset() // 清空缓冲，下次读取会重新请求
	return u.FilePos, nil
}

// fillBuffer 从远程文件下载数据
func (u *URL) fillBuffer() error {
	if u.client == nil {
		return errors.New("http client not initialized")
	}
	req, err := http.NewRequest("GET", u.url, nil)
	if err != nil {
		return err
	}
	// 支持 Range 请求
	rangeHeader := "bytes=" + strconv.FormatInt(u.FilePos, 10) + "-" + strconv.FormatInt(u.FilePos+65535, 10)
	req.Header.Set("Range", rangeHeader)

	resp, err := u.client.Do(req)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	buf := new(bytes.Buffer)
	n, err := buf.ReadFrom(resp.Body)
	if err != nil && err != io.EOF {
		return err
	}
	u.FilePos += n
	u.buf = buf
	return nil
}

