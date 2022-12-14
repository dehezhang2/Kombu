#include <kombu/volume.h>
#include <QFile>
#include <QDataStream>
#include <filesystem/resolver.h>
#include <nori/shape.h>

#if defined(PLATFORM_LINUX) || defined(PLATFORM_MACOS)
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#if defined(PLATFORM_WINDOWS)
#include <windows.h>
#endif

NORI_NAMESPACE_BEGIN
class GridVolume : public Volume{
public:
    GridVolume(const PropertyList &props){
        m_filename = getFileResolver()->resolve(props.getString("filename")).str();

		QString m_filename_q = QString::fromStdString(m_filename);
        QByteArray filename = m_filename_q.toLocal8Bit();
		QFile file(m_filename_q);

		if (!file.exists())
			throw NoriException("The grid file does not exist!");

		/* Parse the file header */
		file.open(QIODevice::ReadOnly);
		QDataStream stream(&file);
		stream.setByteOrder(QDataStream::LittleEndian);

		qint8 header[3], version; qint32 type;
		stream >> header[0] >> header[1] >> header[2] >> version >> type;

		if (memcmp(header, "VOL", 3) != 0 || version != 3)
			throw NoriException("This is not a valid volume data file!");

		stream >> m_res.x() >> m_res.y() >> m_res.z();
		file.close();


		cout << "Mapping \"" << filename.data() << "\" (" << m_res.x()
			<< "x" << m_res.y() << "x" << m_res.z() << ") into memory .." << endl;
		m_stepSize = -1;
		m_fileSize = (size_t) file.size();
		#if defined(PLATFORM_LINUX) || defined(PLATFORM_MACOS)
			int fd = open(filename.data(), O_RDONLY);
			if (fd == -1)
				throw NoriException("Could not open grid file!");
			m_data = (float *) mmap(NULL, m_fileSize, PROT_READ, MAP_SHARED, fd, 0);
			if (m_data == NULL)
				throw NoriException("mmap(): failed.");
			if (close(fd) != 0)
				throw NoriException("close(): unable to close file descriptor!");
		#elif defined(PLATFORM_WINDOWS)
			m_file = CreateFileA(filename.data(), GENERIC_READ, 
				FILE_SHARE_READ, NULL, OPEN_EXISTING, 
				FILE_ATTRIBUTE_NORMAL, NULL);
			if (m_file == INVALID_HANDLE_VALUE)
				throw NoriException(QString("Could not open \"%1\"!").arg(m_filename_q));
			m_fileMapping = CreateFileMapping(m_file, NULL, PAGE_READONLY, 0, 0, NULL);
			if (m_fileMapping == NULL)
				throw NoriException("CreateFileMapping(): failed.");
			m_data = (float *) MapViewOfFile(m_fileMapping, FILE_MAP_READ, 0, 0, 0);
			if (m_data == NULL)
				throw NoriException("MapViewOfFile(): failed.");
		#endif

		m_data += 12; // Shift past the header
	}

    virtual ~GridVolume(){
        if (m_data) {
			m_data -= 12;

			cout << "Unmapping \"" << m_filename << "\" from memory.." << endl;
			#if defined(PLATFORM_LINUX) || defined(PLATFORM_MACOS)
				int retval = munmap(m_data, m_fileSize);
				if (retval != 0)
					throw NoriException("munmap(): unable to unmap memory!");
			#elif defined(PLATFORM_WINDOWS)
				if (!UnmapViewOfFile(m_data))
					throw NoriException("UnmapViewOfFile(): unable to unmap memory region");
				if (!CloseHandle(m_fileMapping))
					throw NoriException("CloseHandle(): unable to close file mapping!");
				if (!CloseHandle(m_file))
					throw NoriException("CloseHandle(): unable to close file");
			#endif
		}
    }

    float lookup(Shape* shape, const Point3f &_p) {
		// Vector3f extents = shape->getBoundingBox().getExtents(), min = shape->getBoundingBox().min;
		Vector3f extents = shape->getLocalExtent(), min = shape->getLocalMin();
		Point3f p = _p - min;
		p = Point3f(p[0]/extents[0] * (m_res[0]-1), p[1]/extents[1] * (m_res[1]-1), p[2]/extents[2] * (m_res[2]-1));
		const int x1 = (int)std::floor(p.x()),
		      y1 = (int)std::floor(p.y()),
		      z1 = (int)std::floor(p.z()),
		      x2 = x1+1, y2 = y1+1, z2 = z1+1;
		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x() ||
		    y2 >= m_res.y() || z2 >= m_res.z())
		    return 0;
		const float fx = p.x() - x1, fy = p.y() - y1, fz = p.z() - z1,
		        _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;
		const float *floatData = (float *) m_data;
		const float
		    d000 = floatData[(z1*m_res.y() + y1)*m_res.x() + x1],
		    d001 = floatData[(z1*m_res.y() + y1)*m_res.x() + x2],
		    d010 = floatData[(z1*m_res.y() + y2)*m_res.x() + x1],
		    d011 = floatData[(z1*m_res.y() + y2)*m_res.x() + x2],
		    d100 = floatData[(z2*m_res.y() + y1)*m_res.x() + x1],
		    d101 = floatData[(z2*m_res.y() + y1)*m_res.x() + x2],
		    d110 = floatData[(z2*m_res.y() + y2)*m_res.x() + x1],
		    d111 = floatData[(z2*m_res.y() + y2)*m_res.x() + x2];
		float ret = ((d000*_fx + d001*fx)*_fy +
                (d010*_fx + d011*fx)*fy)*_fz +
               ((d100*_fx + d101*fx)*_fy +
                (d110*_fx + d111*fx)*fy)*fz;
		return ret;
    }

    float getStepSize(Shape* shape){
		if(m_stepSize == -1){
			Vector3f extents = shape->getBoundingBox().getExtents();
			m_stepSize = std::numeric_limits<float>::infinity();
			for (int i=0; i<3; ++i)
				m_stepSize = std::min(m_stepSize, 0.5f * extents[i] / (float) (m_res[i]-1));
		}
		
        return m_stepSize;
    }

    float getMaximumValue() const {
        return 1.0f;
    }

    std::string toString() const {
        return tfm::format(
                    "GridVolume[\n"
                    "  res = %s,\n"
                    "]",
                    m_res);
    }


protected:
    std::string m_filename;
	size_t m_fileSize;
	float *m_data;
    Vector3i m_res;
	float m_stepSize;
};
NORI_REGISTER_CLASS(GridVolume, "gridvolume");
NORI_NAMESPACE_END