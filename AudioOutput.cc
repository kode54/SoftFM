/*
 *  Audio output handling for SoftFM
 *
 *  Copyright (C) 2013, Joris van Rantwijk.
 *
 *  .WAV file writing by Sidney Cadot,
 *  adapted for SoftFM by Joris van Rantwijk.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, see http://www.gnu.org/licenses/gpl-2.0.html
 */

#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <fcntl.h>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <algorithm>

#ifndef __APPLE__
#include <alsa/asoundlib.h>
#endif

#include "SoftFM.h"
#include "AudioOutput.h"

using namespace std;


/* ****************  class AudioOutput  **************** */

// Encode a list of samples as signed 16-bit little-endian integers.
void AudioOutput::samplesToInt16(const SampleVector& samples,
                                 vector<uint8_t>& bytes)
{
    bytes.resize(2 * samples.size());

    SampleVector::const_iterator i = samples.begin();
    SampleVector::const_iterator n = samples.end();
    vector<uint8_t>::iterator k = bytes.begin();

    while (i != n) {
        Sample s = *(i++);
        s = max(Sample(-1.0), min(Sample(1.0), s));
        long v = lrint(s * 32767);
        unsigned long u = v;
        *(k++) = u & 0xff;
        *(k++) = (u >> 8) & 0xff;
    }
}


/* ****************  class RawAudioOutput  **************** */

// Construct raw audio writer.
RawAudioOutput::RawAudioOutput(const string& filename)
{
    if (filename == "-") {

        m_fd = STDOUT_FILENO;

    } else {

        m_fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
        if (m_fd < 0) {
            m_error  = "can not open '" + filename + "' (" +
                       strerror(errno) + ")";
            m_zombie = true;
            return;
        }

    }
}


// Destructor.
RawAudioOutput::~RawAudioOutput()
{
    // Close file descriptor.
    if (m_fd >= 0 && m_fd != STDOUT_FILENO) {
        close(m_fd);
    }
}


// Write audio data.
bool RawAudioOutput::write(const SampleVector& samples)
{
    if (m_fd < 0)
        return false;

    // Convert samples to bytes.
    samplesToInt16(samples, m_bytebuf);

    // Write data.
    size_t p = 0;
    size_t n = m_bytebuf.size();
    while (p < n) {

        ssize_t k = ::write(m_fd, m_bytebuf.data() + p, n - p);
        if (k <= 0) {
            if (k == 0 || errno != EINTR) {
                m_error = "write failed (";
                m_error += strerror(errno);
                m_error += ")";
                return false;
            }
        } else {
            p += k;
        }
    }

    return true;
}


/* ****************  class WavAudioOutput  **************** */

// Construct .WAV writer.
WavAudioOutput::WavAudioOutput(const std::string& filename,
                               unsigned int samplerate,
                               bool stereo)
  : numberOfChannels(stereo ? 2 : 1)
  , sampleRate(samplerate)
{
    m_stream = fopen(filename.c_str(), "wb");
    if (m_stream == NULL) {
        m_error  = "can not open '" + filename + "' (" +
                   strerror(errno) + ")";
        m_zombie = true;
        return;
    }

    // Write initial header with a dummy sample count.
    // This will be replaced with the actual header once the WavFile is closed.
    if (!write_header(0x7fff0000)) {
        m_error = "can not write to '" + filename + "' (" +
                  strerror(errno) + ")";
        m_zombie = true;
    }
}


// Destructor.
WavAudioOutput::~WavAudioOutput()
{
    // We need to go back and fill in the header ...

    if (!m_zombie) {

        const unsigned bytesPerSample = 2;

        const long currentPosition = ftell(m_stream);

        assert((currentPosition - 44) % bytesPerSample == 0);

        const unsigned totalNumberOfSamples = (currentPosition - 44) / bytesPerSample;

        assert(totalNumberOfSamples % numberOfChannels == 0);

        // Put header in front

        if (fseek(m_stream, 0, SEEK_SET) == 0) {
            write_header(totalNumberOfSamples);
        }
    }

    // Done writing the file

    if (m_stream) {
        fclose(m_stream);
    }
}


// Write audio data.
bool WavAudioOutput::write(const SampleVector& samples)
{
    if (m_zombie)
        return false;

    // Convert samples to bytes.
    samplesToInt16(samples, m_bytebuf);

    // Write samples to file.
    size_t k = fwrite(m_bytebuf.data(), 1, m_bytebuf.size(), m_stream);
    if (k != m_bytebuf.size()) {
        m_error = "write failed (";
        m_error += strerror(errno);
        m_error += ")";
        return false;
    }

    return true;
}


// (Re)write .WAV header.
bool WavAudioOutput::write_header(unsigned int nsamples)
{
    const unsigned bytesPerSample = 2;
    const unsigned bitsPerSample  = 16;

    enum wFormatTagId
    {
        WAVE_FORMAT_PCM        = 0x0001,
        WAVE_FORMAT_IEEE_FLOAT = 0x0003
    };

    assert(nsamples % numberOfChannels == 0);

    // synthesize header

    uint8_t wavHeader[44];

    encode_chunk_id    (wavHeader +  0, "RIFF");
    set_value<uint32_t>(wavHeader +  4, 36 + nsamples * bytesPerSample);
    encode_chunk_id    (wavHeader +  8, "WAVE");
    encode_chunk_id    (wavHeader + 12, "fmt ");
    set_value<uint32_t>(wavHeader + 16, 16);
    set_value<uint16_t>(wavHeader + 20, WAVE_FORMAT_PCM);
    set_value<uint16_t>(wavHeader + 22, numberOfChannels);
    set_value<uint32_t>(wavHeader + 24, sampleRate                                    ); // sample rate
    set_value<uint32_t>(wavHeader + 28, sampleRate * numberOfChannels * bytesPerSample); // byte rate
    set_value<uint16_t>(wavHeader + 32,              numberOfChannels * bytesPerSample); // block size
    set_value<uint16_t>(wavHeader + 34, bitsPerSample);
    encode_chunk_id    (wavHeader + 36, "data");
    set_value<uint32_t>(wavHeader + 40, nsamples * bytesPerSample);

    return fwrite(wavHeader, 1, 44, m_stream) == 44;
}


void WavAudioOutput::encode_chunk_id(uint8_t * ptr, const char * chunkname)
{
    for (unsigned i = 0; i < 4; ++i)
    {
        assert(chunkname[i] != '\0');
        ptr[i] = chunkname[i];
    }
    assert(chunkname[4] == '\0');
}


template <typename T>
void WavAudioOutput::set_value(uint8_t * ptr, T value)
{
    for (size_t i = 0; i < sizeof(T); ++i)
    {
        ptr[i] = value & 0xff;
        value >>= 8;
    }
}


#ifndef __APPLE__
/* ****************  class AlsaAudioOutput  **************** */

// Construct ALSA output stream.
AlsaAudioOutput::AlsaAudioOutput(const std::string& devname,
                                 unsigned int samplerate,
                                 bool stereo)
{
    m_pcm = NULL;
    m_nchannels = stereo ? 2 : 1;

    int r = snd_pcm_open(&m_pcm, devname.c_str(),
                         SND_PCM_STREAM_PLAYBACK, SND_PCM_NONBLOCK);

    if (r < 0) {
        m_error = "can not open PCM device '" + devname + "' (" +
                  strerror(-r) + ")";
        m_zombie = true;
        return;
    }

    snd_pcm_nonblock(m_pcm, 0);

    r = snd_pcm_set_params(m_pcm,
                           SND_PCM_FORMAT_S16_LE,
                           SND_PCM_ACCESS_RW_INTERLEAVED,
                           m_nchannels,
                           samplerate,
                           1,               // allow soft resampling
                           500000);         // latency in us

    if (r < 0) {
        m_error = "can not set PCM parameters (";
        m_error += strerror(-r);
        m_error += ")";
        m_zombie = true;
    }
}


// Destructor.
AlsaAudioOutput::~AlsaAudioOutput()
{
    // Close device.
    if (m_pcm != NULL) {
        snd_pcm_close(m_pcm);
    }
}


// Write audio data.
bool AlsaAudioOutput::write(const SampleVector& samples)
{
    if (m_zombie)
        return false;

    // Convert samples to bytes.
    samplesToInt16(samples, m_bytebuf);

    // Write data.
    unsigned int p = 0;
    unsigned int n = samples.size() / m_nchannels;
    unsigned int framesize = 2 * m_nchannels;
    while (p < n) {

        int k = snd_pcm_writei(m_pcm,
                               m_bytebuf.data() + p * framesize, n - p);
        if (k < 0) {
            m_error = "write failed (";
            m_error += strerror(errno);
            m_error += ")";
            // After an underrun, ALSA keeps returning error codes until we
            // explicitly fix the stream.
            snd_pcm_recover(m_pcm, k, 0);
            return false;
        } else {
            p += k;
        }
    }

    return true;
}
#endif

#ifdef __APPLE__
static const uint32_t DEFAULT_CHUNK_MS = 60;

CoreAudioOutput::CoreAudioOutput(unsigned int sampleRate, bool stereo) :
    m_nchannels(stereo ? 2 : 1), m_audioQueue(NULL), m_sampleRate(sampleRate)
{
    const unsigned int bufferSize = 2048;
    const unsigned int audioLatencyFrames = m_sampleRate * DEFAULT_CHUNK_MS / 1000;
    m_bufferByteSize = bufferSize << m_nchannels;
    // Number of buffers should be ceil(audioLatencyFrames / bufferSize)
    m_numberOfBuffers = (audioLatencyFrames + bufferSize - 1) / bufferSize;
    m_buffers = new AudioQueueBufferRef[m_numberOfBuffers];

    pthread_mutex_init( &m_mutex, NULL );

    AudioStreamBasicDescription dataFormat = {static_cast<Float64>(m_sampleRate), kAudioFormatLinearPCM, kAudioFormatFlagIsSignedInteger | kAudioFormatFlagsNativeEndian, 2 * m_nchannels, 1, 2 * m_nchannels, m_nchannels, 16, 0};
    OSStatus res = AudioQueueNewOutput(&dataFormat, renderOutputBuffer, this, NULL, NULL, 0, &m_audioQueue);
    if (res || m_audioQueue == NULL) {
        return;
    }

    for (uint i = 0; i < m_numberOfBuffers; i++) {
        res = AudioQueueAllocateBuffer(m_audioQueue, m_bufferByteSize, m_buffers + i);
        if (res || m_buffers[i] == NULL) {
            res = AudioQueueDispose(m_audioQueue, true);
            m_audioQueue = NULL;
            return;
        }
        m_buffers[i]->mAudioDataByteSize = m_bufferByteSize;
        // Prime the buffer allocated
        renderOutputBuffer(this, NULL, m_buffers[i]);
    }

    UInt32 numFramesPrepared;
    AudioQueuePrime(m_audioQueue, 0, &numFramesPrepared);

    res = AudioQueueStart(m_audioQueue, NULL);
    if (res) {
        res = AudioQueueDispose(m_audioQueue, true);
        m_audioQueue = NULL;
    }
}

void CoreAudioOutput::renderOutputBuffer(void *userData, AudioQueueRef queue, AudioQueueBufferRef buffer) {
    CoreAudioOutput *stream = (CoreAudioOutput *)userData;
    if (queue == NULL) {
        // Priming the buffers, skip timestamp handling
        queue = stream->m_audioQueue;
    }

    uint byteCount = buffer->mAudioDataByteSize;

    pthread_mutex_lock( &stream->m_mutex );
    if (byteCount > stream->m_bytebuf.size()) {
        memset(buffer->mAudioData, 0, byteCount);
    }
    else {
        memcpy( buffer->mAudioData, &stream->m_bytebuf[0], byteCount );
        stream->m_bytebuf.erase(stream->m_bytebuf.begin(), stream->m_bytebuf.begin() + byteCount);
    }
    pthread_mutex_unlock( &stream->m_mutex );

    AudioQueueEnqueueBuffer(queue, buffer, 0, NULL);
}

CoreAudioOutput::~CoreAudioOutput() {
    if (m_audioQueue) {
        AudioQueueDispose(m_audioQueue, true);
    }
}

bool CoreAudioOutput::write(const SampleVector& samples) {
    // Convert samples to bytes.
    samplesToInt16(samples, m_convbuf);

    pthread_mutex_lock( &m_mutex );
    while ( m_bytebuf.size() >= 8192 ) {
        pthread_mutex_unlock( &m_mutex );
        usleep( 10000 );
        pthread_mutex_lock( &m_mutex );
    }

    m_bytebuf.insert( m_bytebuf.end(), m_convbuf.begin(), m_convbuf.end() );
    pthread_mutex_unlock( &m_mutex );

    return true;
}

#endif

/* end */
