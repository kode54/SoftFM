
#include <cstring>
#include <rtl-sdr.h>

#include "RtlSdrSource.h"

using namespace std;


// Open RTL-SDR device.
RtlSdrSource::RtlSdrSource(int dev_index)
    : m_dev(0)
    , m_block_length(default_block_length)
{
    int r;

    r = rtlsdr_open(&m_dev, dev_index);
    if (r < 0) {
        m_error =  "Failed to open RTL-SDR device (";
        m_error += strerror(-r);
        m_error += ")";
    }
}


// Close RTL-SDR device.
RtlSdrSource::~RtlSdrSource()
{
    if (m_dev)
        rtlsdr_close(m_dev);
}


// Configure RTL-SDR tuner and prepare for streaming.
bool RtlSdrSource::configure(uint32_t sample_rate,
                             uint32_t frequency,
                             int tuner_gain,
                             int block_length)
{
    int r;

    if (!m_dev)
        return false;

    r = rtlsdr_set_sample_rate(m_dev, sample_rate);
    if (r < 0) {
        m_error = "rtlsdr_set_sample_rate failed";
        return false;
    }

    r = rtlsdr_set_center_freq(m_dev, frequency);
    if (r < 0) {
        m_error = "rtlsdr_set_center_freq failed";
        return false;
    }

    if (tuner_gain < 0) {
        r = rtlsdr_set_tuner_gain_mode(m_dev, 0);
        if (r < 0) {
            m_error = "rtlsdr_set_tuner_gain_mode could not set automatic gain";
            return false;
        }
    } else {
        r = rtlsdr_set_tuner_gain_mode(m_dev, 1);
        if (r < 0) {
            m_error = "rtlsdr_set_tuner_gain_mode could not set manual gain";
            return false;
        }

        r = rtlsdr_set_tuner_gain(m_dev, tuner_gain);
        if (r < 0) {
            m_error = "rtlsdr_set_tuner_gain failed";
            return false;
        }
    }

    // set block length
    m_block_length = (block_length < 4096) ? 4096 :
                     (block_length > 1024 * 1024) ? 1024 * 1024 :
                     block_length;
    m_block_length -= m_block_length % 4096;

    // reset buffer to start streaming
    if (rtlsdr_reset_buffer(m_dev) < 0) {
        m_error = "rtlsdr_reset_buffer failed";
        return false;
    }

    return true;
}


// Return current sample frequency in Hz.
uint32_t RtlSdrSource::get_sample_rate()
{
    return rtlsdr_get_sample_rate(m_dev);
}


// Return current center frequency in Hz.
uint32_t RtlSdrSource::get_frequency()
{
    return rtlsdr_get_center_freq(m_dev);
}


// Return current tuner gain in dB.
double RtlSdrSource::get_tuner_gain()
{
    return 0.1 * rtlsdr_get_tuner_gain(m_dev);
}


// Return a list of supported tuner gain settings in dB.
vector<double> RtlSdrSource::get_tuner_gains()
{
    vector<double> result;

    int num_gains = rtlsdr_get_tuner_gains(m_dev, NULL);
    if (num_gains <= 0)
        return result;

    int gains[num_gains];
    if (rtlsdr_get_tuner_gains(m_dev, gains) != num_gains)
        return result;

    result.reserve(num_gains);
    for (int i = 0; i < num_gains; i++)
        result.push_back(0.1 * gains[i]);

    return result;
}


// Fetch a bunch of samples from the device.
bool RtlSdrSource::get_samples(IQSampleVector& samples)
{
    int r, n_read;

    if (!m_dev)
        return false;

    vector<uint8_t> buf(2 * m_block_length);

    r = rtlsdr_read_sync(m_dev, buf.data(), 2 * m_block_length, &n_read);
    if (r < 0) {
        m_error = "rtlsdr_read_sync failed";
        return false;
    }

    if (n_read != 2 * m_block_length) {
        m_error = "short read, samples lost";
        return false;
    }

    samples.resize(m_block_length);
    for (int i = 0; i < m_block_length; i++) {
        int32_t re = buf[2*i];
        int32_t im = buf[2*i+1];
        samples[i] = IQSample( (re - 128) / IQSample::value_type(128),
                               (im - 128) / IQSample::value_type(128) );
    }

    return true;
}


// Return a list of supported devices.
vector<string> RtlSdrSource::get_device_names()
{
    vector<string> result;

    int device_count = rtlsdr_get_device_count();
    if (device_count <= 0)
        return result;

    result.reserve(device_count);
    for (int i = 0; i < device_count; i++) {
        result.push_back(string(rtlsdr_get_device_name(i)));
    }

    return result;
}

/* end */